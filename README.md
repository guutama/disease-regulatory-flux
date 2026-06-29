# Disease-regulatory flux map

A Nextflow pipeline that builds a **disease-regulatory flux map**: it splits each gene's
disease association into an own (cis) part and a network (trans) part, and attributes the
trans part to the specific upstream regulators that deliver it.

The pipeline is **data-agnostic** — it works on any cohort by changing configuration, not
code. This README describes the inputs you provide; pipeline stages are documented as they
are added.

---

## Inputs

The pipeline runs **per tissue** and needs three matched reference files. All are
tab-separated tables and may be gzipped (`.gz`). Column names and the delimiter are
configurable; the defaults are shown below.

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

A variant-by-sample matrix of **alternate-allele dosages in [0, 2]**. The leading
column(s) identify the variant (see *Variant identifiers* below); the remaining columns are
the samples — the **same samples** as in the expression matrix.

```
rsID      S1     S2     S3
rs1001    0.0    1.0    2.0
rs1002    1.0    0.0    1.0
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

### Variant identifiers

The genotype and eQTL files must identify variants the **same way**:

- if **both** files have an `rsID` column, variants are matched by rsID;
- otherwise **both** must have `chromosome`, `position`, `ref`, `alt` columns, and variants
  are matched by `chr:pos:ref:alt`.

### How the files connect

| key | links |
|-----|-------|
| sample columns | expression ↔ genotype |
| `gene_id`      | expression ↔ eQTL |
| variant id     | genotype ↔ eQTL |

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

Then run:

```
nextflow run flux.nf --samplesheet samplesheet.csv --outdir results
```

Outputs are written per tissue under `results/harmonised/`. The column names and field
delimiter are set in `nextflow.config` (or overridden with `--<param>`), so a new dataset
needs no code changes.

---

## First step: harmonisation

The first stage aligns the inputs so the rest of the pipeline can rely on them
(`bin/harmonise_inputs.py`):

- keeps **all** expression genes (every gene is a potential network target);
- reduces the **genotype to the variants present in the eQTL map** (the only variants the
  method uses);
- restricts both expression and genotype to the **samples they share**, in one consistent
  order;
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
inputs using [BioFindr](https://github.com/tmichoel/BioFindr) (`bin/reconstruct_grn.jl`). It
runs the canonical three-function pipeline---`findr` (pairwise causal posteriors, each edge
oriented by the regulator's instrument), `globalfdr!` (keep edges at the chosen
false-discovery rate) and `dagfindr!` (greedy cycle removal)---to produce a directed acyclic
network in which every gene is a candidate target and every edge starts from an
instrument-anchored regulator. Per tissue it writes:

- `<tissue>.grn_edges.csv` — the network edges (`Source`, `Target`, edge probability and
  q-value);
- `<tissue>.grn_edges_all.csv.gz` — every inferred edge with its q-value, before the FDR cut,
  so a network at a different FDR can be derived without recomputing the posteriors;
- `<tissue>.grn_qc.tsv` — node, edge, regulator, target and degree counts;

under `results/grn/`. This step uses Julia and BioFindr (pinned in `julia_project/`); the
test combination, FDR and mixture-fit method are set by `grn_combination`, `grn_fdr` and
`grn_findr_method` in `nextflow.config`.

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
