# Analysis code

Scripts that produce the figures and tables in the disease-regulatory-flux manuscript, kept
separate from the reproducible pipeline (`bin/`, `modules/`, `flux.nf`). Each script is the code
as run, with a module docstring saying what it makes, which inputs it reads, and how to run it.

```
analysis/
  figures/      one script per manuscript figure (main text and supplementary)
  tables/       builder for the release supplementary tables (S1-S9)
  validation/   apply the published predictors to an external GWAS or cohort
```

## Data availability

STARNET is individual-level, controlled-access data (dbGaP/EGA); the raw genotypes, expression
matrices and cohort phenotypes require data-access approval and are **not** included in this
repository. The scripts read the pipeline's per-gene, per-edge and per-tissue *summary* result
files, and write figures — no individual-level data is committed or emitted. The one script that
needs the restricted cohort workbook (`plot_covariate_check.py`) takes its path as a required
argument and has no default, so it never hard-codes a controlled-access location.

## Environment variables

The scripts resolve their inputs and outputs through three variables, all with repo-local
defaults so they run from a checkout without configuration:

| Variable | Purpose | Default |
|----------|---------|---------|
| `FLUX_RESULTS` | directory of pipeline result files (the inputs) | `results_cv/` |
| `FLUX_FIGURES` | directory the figures are written to | `results/figures/` |
| `GENE_ANNOT` | GENCODE v19 gene annotation TSV, for HGNC symbol labels | — (required by the scripts that label genes) |

```
GENE_ANNOT=/path/to/gencode.v19.genes.tsv \
FLUX_RESULTS=/path/to/results \
python analysis/figures/plot_main_readouts.py
```

Every figure lands in one shared directory (`FLUX_FIGURES`, default `results/figures/`), so a full
figure set is produced by running the scripts in any order. Figures are Arial, 300 DPI.

## Main-text figures

| Script | Figure | What it shows |
|--------|--------|---------------|
| `figures/plot_alltargets_summary.py` | `fig_alltargets_summary` | transcriptome-wide expression-prediction summary |
| `figures/plot_main_readouts.py` | `fig_main_readouts` | six-panel disease-role readouts (needs `GENE_ANNOT`) |
| `figures/plot_flux_map_grn.py` | `fig_flux_map_grn` | disease-regulatory flux map with the full GRN cascade (needs `GENE_ANNOT`) |
| `figures/plot_flux_validation.py` | `fig_flux_validation` | four-panel flux-map validation (OmniPath, negative control) |
| `figures/plot_functional_axis.py` | `fig_functional_axis` | cardiovascular functional-relevance enrichment (needs `GENE_ANNOT`) |
| `figures/plot_assoc_combined.py` | `fig_assoc_combined` | combined CAD transcriptome-wide association (Manhattan via `utils/`, needs `GENE_ANNOT`) |

## Supplementary figures

These reproduce the figures in `paper/flux_map_supplementary.tex`.

| Script | Figure | What it shows |
|--------|--------|---------------|
| `figures/plot_covariate_check.py` | `fig_covariate_check`, `fig_covariate_pc_scatter` | expression covariate-QC: PC-vs-covariate heatmap and PC scatter (needs `--pheno-xls`) |
| `figures/plot_coexpression_density.py` | `fig_coexpression_density` | background co-expression versus reconstructed-network density per tissue |
| `figures/plot_grn_vs_classical_5methods.py` | `fig_grn_vs_classical_5methods` | five-method prediction comparison (three GRN priors vs two classical cis-only) |
| `figures/compare_traingrn_leakage.py` | `fig_leakage_traingrn` | data-leakage control: full-cohort versus 80%-train re-run |
| `figures/plot_flux_omnipath.py` | `fig_flux_omnipath_network` | OmniPath corroboration of the flux-map regulator→disease relations |

## Supplementary tables

`tables/build_supp_tables.py` assembles the release supplementary tables (S1-S9) from the
per-tissue result files — expression-model summaries (S1), CAD association (S2), GRN edges (S3),
flux edges (S4), flux node roles (S5), FinnGen replication (S6), cardiovascular gene-set
membership (S7), functional enrichment (S8), and per-SNP predictor weights with dbSNP alleles and
regulator provenance (S9) — plus a README data dictionary. Every table carries gene-, edge- and
tissue-level summaries only. S9 (with alleles) is built when `--variant-map` is supplied.

```
GENE_ANNOT=/path/to/gencode.v19.genes.tsv \
python analysis/tables/build_supp_tables.py \
    --outdir supplementary_tables \
    --variant-map <variant annotation TSV> --variant-map-id-col original_id
```

## External validation

Both validators reuse the published predictor weights (Supplementary Table S9) with no retraining
and no access to the individual-level STARNET data.

**Summary-statistic TWAS on a new trait** — `validation/twas_from_gwas.py` applies the weights to
any GWAS to obtain a per-gene, per-tissue TWAS statistic, decomposed into cis and trans
components. Align the GWAS to the S9 ALT alleles first (with `bin/harmonise_gwas.py`), then:

```
python analysis/validation/twas_from_gwas.py \
    --weights   Supplementary_Table_S9_predictor_weights.csv.gz \
    --sigma-g   Supplementary_Table_S2_twas_association.csv \
    --gwas      <trait>.aligned.tsv.gz \
    --out       <trait>.external_twas.csv
```

**Predictor portability in a new cohort** — `validation/validate_predictors_external.py` tests our
LOO-R2 (Supplementary Table S1) in an independent cohort with genotypes and measured expression
from a matching tissue. It predicts expression from dosage (`yhat = sum_j weight_j * dosage_j`,
cis and trans channels together) and reports the external squared correlation against the measured
expression, next to our `loo_r2`. No model is refit.

```
python analysis/validation/validate_predictors_external.py \
    --weights     Supplementary_Table_S9_predictor_weights.csv.gz \
    --loo         Supplementary_Table_S1_expression_models.csv \
    --dosage      <external dosage matrix.tsv.gz> \
    --expression  <external expression matrix.tsv.gz> \
    --tissue      AOR --out <cohort>.external_r2.csv --fig
```

## Shared helpers

Figure scripts import shared helpers from `utils/` (a repo-level package): `gene_labels.py`
(ENSG→HGNC symbols), `figure_data.py` (shared data layer and style), and the Manhattan renderer
(`plot_manhattan.py` + `manhattan_ggplot.R` + `rscript_wrapper.sh`).
