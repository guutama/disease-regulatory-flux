# Disease-regulatory flux map

[![DOI](https://zenodo.org/badge/1282993883.svg)](https://zenodo.org/badge/latestdoi/1282993883)

A Nextflow pipeline that builds a **disease-regulatory flux map**: it splits each gene's
disease association into an own (cis) part and a network (trans) part, and attributes the
trans part to the specific upstream regulators that deliver it.

The pipeline is **data-agnostic** — it works on any cohort by changing configuration, not
code. This README describes the inputs you provide; pipeline stages are documented as they
are added.

---

## Inputs

To run end to end — from reference data to the disease-regulatory flux map — the pipeline needs
**four** inputs: three matched **per-tissue** reference files (expression, genotype and cis-eQTL),
and one or more **GWAS summary-statistics** files for the trait(s) of interest. All are
tab-separated tables and may be gzipped (`.gz`). Column names and the delimiter are configurable;
the defaults are shown below.

### 1. Expression matrix — `expression.tsv[.gz]`

A gene-by-sample matrix of **analysis-ready** expression (already normalised and
covariate-adjusted — the pipeline uses the values as given). The first column is the gene
id; the remaining columns are the samples.

```
gene_id     S1       S2       S3
ENSG0001    0.42    -1.10     0.05
ENSG0002   -0.88     0.31     1.22
ENSG0003    0.10     0.07    -0.44
```

All genes are kept: every expressed gene can be a **target** in the regulatory network.

### 2. Genotype matrix — `genotype.tsv[.gz]`

A variant-by-sample matrix of **alternate-allele dosages in [0, 2]**. The leading columns
identify the variant and **must include `chromosome`, `position`, `ref` and `alt`** (the
alternate allele the dosage counts), optionally with an `rsID`; the remaining columns are
the samples — the **same samples** as in the expression matrix. The allele columns are
carried through the whole pipeline so a GWAS can later be aligned to the alternate allele
(see *Variant identifiers* below).

```
rsID      chromosome  position  ref  alt   S1     S2     S3
rs1001    1           1001      A    G     0.0    1.0    2.0
rs1002    1           1002      C    T     1.0    0.0    1.0
```

### 3. eQTL map — `eqtl.tsv[.gz]`

A long table of **cis-eQTL associations** linking variants to genes — one row per
gene–variant pair, with effect size, standard error and p-value. This is what marks a gene
as a possible **regulator** (only genes with a cis-eQTL can regulate others).

```
gene_id     rsID      beta     se      pvalue
ENSG0001    rs1001    0.31    0.06    2.1e-7
ENSG0001    rs1002   -0.12    0.05    3.0e-3
ENSG0002    rs1002    0.27    0.06    8.0e-6
```

### 4. GWAS summary statistics — `gwas.tsv[.gz]`

Genome-wide association summary statistics for the trait — one row per variant, with the effect
allele, the other allele, effect size, standard error and p-value. This is the disease signal the
pipeline carries through the network to build the flux map. Unlike the three files above it is
**not** per tissue: provide **one file per trait**, and the same GWAS is used across all tissues.

```
chromosome  position  effect_allele  other_allele  beta     se      pvalue    n
1           1001       G              A             0.031    0.006   2.0e-7    500000
1           1002       T              C            -0.012    0.005   3.0e-3    500000
```

Variants are matched to the pipeline by `chromosome` and `position` and re-aligned to the
genotype's alternate allele (Step 8). The column names default to the GWAS-Catalog harmonised
layout and can be overridden per cohort.

### Variant identifiers

The genotype and eQTL files must identify variants the **same way**:

- if **both** files have an `rsID` column, variants are matched by rsID;
- otherwise **both** must have `chromosome`, `position`, `ref`, `alt` columns, and variants
  are matched by `chr:pos:ref:alt`.

Either way, the **genotype must always carry `chromosome`, `position`, `ref` and `alt`** (even
when matching by rsID). These alleles are required and propagated through every stage so the
GWAS-harmonisation step can align each variant to its alternate allele.

### How the files connect

| key | links |
|-----|-------|
| sample columns | expression ↔ genotype |
| `gene_id`      | expression ↔ eQTL |
| variant id     | genotype ↔ eQTL |
| `chromosome`, `position` | genotype ↔ GWAS |

---

## Running the pipeline

Requires Nextflow ≥ 24.04.2 (and Java 11+).

You provide the pipeline with a **samplesheet**: a small CSV file you write yourself that
lists, for each tissue, where its three input files are. The samplesheet has a header row
with four fixed column names and then **one row per tissue**:

```
tissue,expression,genotype,eqtl
AOR,/data/AOR.expr.tsv.gz,/data/AOR.geno.tsv.gz,/data/AOR.eqtl.tsv.gz
LIV,/data/LIV.expr.tsv.gz,/data/LIV.geno.tsv.gz,/data/LIV.eqtl.tsv.gz
```

| column | what to put in it |
|--------|-------------------|
| `tissue` | a short label for the tissue (your choice); used to name that tissue's outputs |
| `expression` | path to the tissue's expression matrix (the file described above) |
| `genotype` | path to the tissue's genotype dosage matrix |
| `eqtl` | path to the tissue's cis-eQTL table |

Add one row for every tissue you want to run. Paths may be absolute or relative to where
you launch the pipeline. The input files themselves are the `.tsv[.gz]` tables described in
*Inputs* above — the samplesheet only points to them.

You also provide a **GWAS samplesheet** — a second CSV with a header and one row per trait,
naming each trait's summary-statistics file:

```
trait,gwas
CAD,/data/cad_gwas.tsv.gz
```

| column | what to put in it |
|--------|-------------------|
| `trait` | a short label for the trait (your choice); used to name that trait's outputs |
| `gwas` | path to the trait's GWAS summary-statistics file (the file described above) |

Then run the whole pipeline end to end:

```
nextflow run flux.nf \
    --samplesheet samplesheet.csv \
    --gwas_samplesheet gwas.csv \
    --outdir results
```

This produces the disease-regulatory flux map per tissue and trait under `results/flux/`. The
column names and field delimiter are set in `nextflow.config` (or overridden with `--<param>`), so
a new dataset needs no code changes. Omitting `--gwas_samplesheet` runs only the expression-model
half of the pipeline and stops after Step 7 (no association or flux map).

---

## First step: harmonisation

The first stage aligns the inputs so the rest of the pipeline can rely on them
(`bin/harmonise_inputs.py`):

- keeps **all** expression genes (every gene is a potential network target);
- reduces the **genotype to the variants present in the eQTL map** (the only variants the
  method uses);
- restricts both expression and genotype to the **samples they share**, in one consistent
  order;
- carries the variant's `chromosome`, `position`, `ref` and `alt` through into the harmonised
  genotype (`variant_id`, `chromosome`, `position`, `ref`, `alt`, then the samples);
- writes a QC report of how many samples, genes and variants were kept.

It stops with a clear error if the files don't share samples, or if no eQTL pair survives
the matching.

---

## Second step: instrument selection

The second stage prepares the inputs for reconstructing the regulatory network
(`bin/select_instruments.py`). A gene can regulate others only if it has a cis-eQTL, and
the network reconstruction uses one genetic instrument per regulator — its **lead** cis-eQTL
(the SNP with the smallest p-value; ties broken by the largest effect size, then the variant
id). From the harmonised trio it writes three per-tissue files:

- `<tissue>.dE.csv` — the **instrument list**: one row per regulator, as a `variant_id,gene_id`
  pair;
- `<tissue>.dG.csv` — the **genotype** of those instruments, samples in rows and one column
  per instrument SNP;
- `<tissue>.dX.csv` — the **expression** of all genes, samples in rows and one column per gene
  (every gene is a candidate target);

plus a QC report of how many samples, genes and instruments were selected. The expression
and genotype files share one sample order, so they line up row by row. The files are written
under `results/grn_inputs/`.

---

## Third step: GRN reconstruction

The third stage reconstructs each tissue's gene regulatory network from the `dX`/`dG`/`dE`
inputs using [findr](https://github.com/lingfeiwang/findr) (`bin/run_findr_py.py`), the
python package around Lingfei Wang's libfindr C core. It runs the canonical three-step
pipeline---pairwise causal posteriors (each edge oriented by the regulator's instrument), a
global false-discovery q-value cut, and greedy cycle removal---to produce a directed acyclic
network in which every gene is a candidate target and every edge starts from an
instrument-anchored regulator. Per tissue it writes:

- `<tissue>.grn_edges.csv` — the network edges (`Source`, `Target`, edge probability and
  q-value);
- `<tissue>.grn_edges_all.csv.gz` — every inferred edge with its q-value, before the FDR cut,
  so a network at a different FDR can be derived without recomputing the posteriors;
- `<tissue>.grn_qc.tsv` — node, edge, regulator, target and degree counts;

under `results/grn/`. The libfindr shared library and its GSL dependency are located via
`findr_libpath` and `findr_gsl_dir`; the test combination and FDR are set by
`grn_combination` and `grn_fdr` in `nextflow.config`.

---

## Fourth step: cis-eQTL feature selection

The fourth stage builds each gene's **cis channel** --- the set of cis-eQTL SNPs used to
predict its expression (`bin/select_cis_features.py`). The harmonised eQTL table is already
cis (the eQTL analysis windowed it), so this step does not re-window: it collects, per gene,
the variants associated with it and drops duplicate pairs. If `cis_max_pvalue` is set in
`nextflow.config` and the table has a p-value column, weaker associations are dropped;
otherwise every association is trusted. Per tissue it writes:

- `<tissue>.cis_features.tsv.gz` — the per-gene cis-SNP set (`gene_id`, `variant_id`);
- `<tissue>.cis_features_summary.tsv` — per gene: number of cis SNPs;
- `<tissue>.cis_features_qc.tsv` — pairs read, kept and dropped;

under `results/cis_features/`.

---

## Fifth step: LD filtering

The fifth stage thins each gene's cis-SNP set by linkage disequilibrium and emits a hard-call
genotype matrix for the SNPs that survive (`bin/ld_prune.py`). Neighbouring SNPs are often
strongly correlated, so the raw cis set is redundant. Per gene, the SNPs are ordered from the
most to the least significant cis-eQTL and walked in turn: a SNP is kept unless its squared
correlation (on the tissue's samples, using 0/1/2 hard calls) with an already-kept SNP
exceeds `ld_r2_threshold` (default `0.7`, set in `nextflow.config`). Constant SNPs are
dropped. Per tissue it writes:

- `<tissue>.cis_snps_pruned.tsv.gz` — the kept cis SNPs per gene (`gene_id`, `variant_id`);
- `<tissue>.genotype_012.tsv.gz` — a 0/1/2 hard-call matrix (`variant_id`, `chromosome`,
  `position`, `ref`, `alt`, then samples) for the kept SNPs;
- `<tissue>.ld_prune_summary.tsv` — per gene: cis SNPs in, kept and dropped;
- `<tissue>.ld_prune_qc.tsv` — tissue-level totals;

under `results/ld_pruned/`. This step uses NumPy.

---

## Sixth step: trans-feature construction

The sixth stage builds each gene's **trans channel** --- the genetics of its upstream
regulators (`bin/trans_features.py`). Starting from each gene it walks the reconstructed
network backwards: to the gene's direct regulators (hop 1), their regulators (hop 2), and so
on up to `trans_max_hop` hops (default `2`, set in `nextflow.config`). For every regulator
reached it collects that regulator's LD-pruned cis SNPs; those SNPs are the gene's trans
features. A regulator reachable at more than one hop is reported at each. Per tissue it writes:

- `<tissue>.trans_features.tsv.gz` — one row per (gene, hop, regulator, variant):
  `gene_id`, `hop`, `source_gene_id`, `variant_id`;
- `<tissue>.trans_features_summary.tsv` — per gene: number of regulators and unique trans SNPs;
- `<tissue>.trans_features_qc.tsv` — tissue-level totals;

under `results/trans_features/`.

---

## Seventh step: expression models

The seventh stage fits a channel-aware Bayesian expression model for every gene
(`bin/fit_expression_models.py`). Each gene's expression is predicted from two genetic
channels --- its own cis SNPs (the cis channel) and its upstream regulators' cis SNPs (the
trans channel) --- and the model is fit in up to three configurations: `cis_only`,
`trans_only` and `cis_trans`. A gene with both channels is fit in all three; a gene with one
channel is fit in that one. Each configuration is fit by MCMC (NumPyro); the configuration with
the highest expected log predictive density (elpd) is selected, and the gene is **predictable**
when the selected model's leave-one-out R² reaches the threshold (`0.01`).

### How predictors are scored: whole-cohort vs held-out (`train_fraction`)

How a configuration is scored is controlled by `train_fraction` in `nextflow.config` (or
`--train_fraction` on the command line):

- **`train_fraction = 1` (default)** — the whole cohort is used for both fitting and scoring, and
  each configuration is scored by PSIS leave-one-out cross-validation (ArviZ). This is the right
  choice for small reference cohorts (a few hundred samples), where a held-out fold would be too
  small to estimate a per-gene R² without dominating noise.
- **`train_fraction < 1`** — an extra step (`SPLIT_SAMPLES`) partitions each tissue's samples into
  a train and a test fold (reproducibly, seeded by `split_seed`). The model is then fit on the
  train fold and scored **out-of-sample** on the held-out test fold: the reported `loo_r2` and
  `elpd` columns become the genuine held-out test-set R² and log predictive density, and the
  written weights are the train-only fit. Use this to check that the predictive accuracy is not
  inflated by fitting and scoring on the same samples. For example, an 80/20 split:

  ```
  nextflow run flux.nf --samplesheet samplesheet.csv --outdir results \
      --train_fraction 0.8 --split_seed 1
  ```

  The per-tissue folds are published under `results/splits/` (`<tissue>.train.samples`,
  `<tissue>.test.samples`), and each tissue's `expr_model_stats.tsv` records `eval_mode`
  (`loo` or `heldout_test`) with the `n_train`/`n_test` counts used.

The prior is set by `expr_model` in `nextflow.config`:

- `bayes_ridge` (default) --- a Gaussian ridge with a separate scale per channel;
- `horseshoe` --- a regularised horseshoe (global-local shrinkage with a slab);
- `bslmm` --- a sparse-plus-polygenic prior (a spike-and-slab of large effects on top of a
  small polygenic background).

The number of MCMC iterations and the seed are set by `expr_mcmc_warmup`, `expr_mcmc_samples`
and `expr_seed`. Per tissue it writes, under `results/expression_models/`:

- `<tissue>.expr_model_metrics.tsv.gz` — per gene: per-config LOO-R² and elpd, the gene class,
  the selected config, its genetic standard deviation and intercept, and whether the gene is
  predictable;
- `<tissue>.expr_model_weights.tsv.gz` — posterior-mean raw-dosage weights (`gene_id`,
  `method`, `config`, `channel`, `variant_id`, `weight`);
- `<tissue>.expr_model_selected.tsv` — the predictable genes and their selected config;
- `<tissue>.expr_model_stats.tsv` — tissue-level counts and mean LOO-R² per channel.

This step requires NumPyro and JAX (MCMC) and ArviZ (PSIS-LOO); see `docs/requirements.md`.

---

## Eighth step: GWAS harmonisation

The eighth stage prepares a GWAS for the association by re-expressing every effect as the
effect of the pipeline's **alternate allele** (`bin/harmonise_gwas.py`). The genotype counts
the alternate allele and the predictor weights are in alternate-dosage units, so the GWAS must
be aligned the same way. For each variant in the pipeline's universe (every tissue's
`genotype_012`, which carries `chromosome`, `position`, `ref` and `alt`), the step finds the
matching GWAS record by chromosome and position and aligns it:

- effect on the alternate allele → `z = +beta/se`; effect on the reference allele → `z = -beta/se`;
- effects reported on the complementary strand are resolved by complementing the GWAS alleles;
- strand-ambiguous SNPs (A/T, C/G) and variants with no GWAS record are dropped.

It writes, per trait, under `results/gwas/`:

- `<trait>.gwas.tsv.gz` — one row per aligned variant: `variant_id`, `z`, `beta`, `se`,
  `p_value`, `n`, with `z` and `beta` on the alternate-allele scale;
- `<trait>.gwas_qc.tsv` — totals (variants aligned, dropped ambiguous, dropped unmatched).

The GWAS column names default to the GWAS-Catalog harmonised layout (`chromosome`, `position`,
`effect_allele`, `other_allele`, `beta`, `se`, `pvalue`, `n`) and can be overridden per cohort.
This step uses only the Python standard library.

This step reads the **GWAS samplesheet** described in *Inputs* and *Running the pipeline*
(`--gwas_samplesheet`). Steps 8–10 are the trait-facing half of the pipeline that turns the GWAS
into the disease-regulatory flux map; when no GWAS samplesheet is given the run stops after Step 7
(expression models only).

---

## Ninth step: TWAS association

The ninth stage tests each gene's genetically predicted expression for association with the
trait (`bin/twas_association.py`). For every gene it combines the predictor's ALT-dosage
weights (Step 7) with the ALT-aligned GWAS (Step 8):

    z_TWAS = ( sum_j  w_j * sd_j * z_gwas_j ) / sigma_g

where `w_j` is SNP `j`'s weight, `sd_j` its genotype standard deviation, `z_gwas_j` the GWAS
z of its alternate allele, and `sigma_g` the spread of the predicted expression over the
genotypes (the LD/variance term — no external panel needed). Restricting the numerator to the
cis or the trans SNPs gives a **cis and a trans component that sum to z_TWAS**, so each gene's
association is split into the part carried by its own cis genetics and the part delivered
through the network. A two-sided p-value comes from the standard normal, and a
Benjamini-Hochberg FDR is applied per tissue.

Per tissue and trait it writes, under `results/association/`:

- `association_<tissue>_<trait>.tsv` — one row per gene: `gene_id`, `config`, SNP counts,
  `z_twas`, `z_cis`, `z_trans`, `sigma_g`, `p_value`, `p_adj`, `tissue`, `trait`.

This step uses NumPy.

---

## Tenth step: disease-regulatory flux map

The tenth stage attributes each disease gene's **trans** signal to the specific upstream
regulators that deliver it (`bin/build_flux_map.py`). Because every trans SNP is a cis-eQTL of
one or more of the gene's GRN ancestors, the gene's trans z decomposes exactly into
per-regulator contributions. For a disease gene `g` (significant at `p_adj < assoc_fdr`) and a
regulator `A` that feeds it:

    flux(A → g) = (1/sigma_g) * sum_{s ∈ E(A) ∩ T(g)}  w_{g,s} * sd_s * z_gwas_s

where `E(A)` is `A`'s cis-eQTL set and `T(g)` is `g`'s trans-feature set. A SNP that instruments
several ancestors **splits its term evenly** among them, so the per-regulator fluxes recover the
gene's trans z exactly: `sum_A flux(A → g) = z_trans(g)`. The flux is **signed** (positive = `A`'s
genetics drive `g` toward the risk-increasing direction) and is an algebraic decomposition, not a
fitted model.

The map keeps the network **paths explicit** rather than collapsing them. Its nodes are the
disease genes plus their hop-1 GRN parents and hop-2 grandparents (the single-hop GRN is recovered
from the hop-1 trans features). Its edges are the true single-hop GRN edges along the paths down to
a disease gene:

- `parent → disease_gene` (hop 1);
- `grandparent → intermediate_parent` (hop 2) — the intermediate parent is kept, not collapsed
  into a `grandparent → disease_gene` shortcut.

Each edge carries the signed flux its source delivers to the disease gene it feeds; a grandparent
reaching a gene through several parents splits its flux evenly across those edges, so the fluxes
toward a gene still sum to `z_trans(g)`. Every parent and grandparent is kept even when it delivers
no flux (its `flux` is left blank).

Per tissue and trait it writes, under `results/flux_map/`:

- `flux_<tissue>_<trait>.tsv` — one row per edge: `source_gene`, `target_gene`, `disease_gene`,
  `hop`, `n_snp`, `flux`, `sign`, `tissue`, `trait`. Together these signed edges are the directed
  disease-regulatory flux network.

This step uses NumPy. Like Steps 8–9, it runs only when a GWAS samplesheet is given.

---

## Dependencies and third-party tools

The pipeline calls a small set of established open-source tools as dependencies; none of their
source is bundled in this repository, so each is fetched from its own distribution and remains
under its own license. The exact versions are pinned in the repository's manifest, which is the
authoritative record: `requirements.txt` for Python.

Pipeline runtime:

- **Nextflow** (>= 24.04.2, Apache-2.0) — workflow orchestration.
- **Python** (>= 3.9) with **NumPy**, **pandas** and **NetworkX**; the instrument-anchored
  causal GRN reconstruction uses **findr** (libfindr, AGPL-3.0) with **GSL**, and the
  Bayesian expression-model step additionally uses **NumPyro**, **JAX** and **ArviZ**.

Analysis and figures (`analysis/`):

- **pandas**, **matplotlib**, **SciPy**, **scikit-learn**, **statsmodels** and **NetworkX**;
  the transcriptome-wide Manhattan panel optionally uses **R** with **ggplot2** and **ggrepel**.

Please cite the tools that ask for it alongside this pipeline — in particular findr
(Wang & Michoel, *PLoS Comput. Biol.* 2017; *Bioinformatics* 2019) and Nextflow
(Di Tommaso et al., *Nat. Biotechnol.* 2017).

## License

This project's own code is released under the MIT License (see [`LICENSE`](LICENSE)).
Third-party dependencies are governed by their respective licenses, noted above.

## Citation

If you use this software, please cite it via the metadata in
[`CITATION.cff`](CITATION.cff) (GitHub's "Cite this repository" button) and the archived
release DOI shown by the badge at the top of this page.
