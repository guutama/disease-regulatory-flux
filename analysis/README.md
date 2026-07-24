# Analysis code

Scripts that produce the figures and derived results in the disease-regulatory-flux manuscript,
kept separate from the reproducible pipeline (`bin/`, `modules/`, `flux.nf`). Each script is the
code as run, with a module docstring saying what it makes, which inputs it reads, and how to run
it.

```
analysis/
  core/           analyses that build the flux network and its cardiovascular enrichment,
                  and fetch the OmniPath reference network
  figure_scripts/ one script per manuscript figure (main text and supplementary)
  validation/     apply the published predictors to an external GWAS or cohort
```

## Data availability

STARNET is individual-level, controlled-access data (dbGaP/EGA); the raw genotypes, expression
matrices and cohort phenotypes require data-access approval and are **not** included in this
repository. The scripts read the pipeline's per-gene, per-edge and per-tissue *summary* result
files and write figures.

## Environment variables

The scripts resolve their inputs and outputs through variables with repo-local defaults so they
run from a checkout without configuration:

| Variable | Purpose | Default |
|----------|---------|---------|
| `FLUX_RESULTS` | directory of pipeline result files (the inputs) | `results_findr_py/` |
| `FLUX_FIGURES` | directory the figures are written to | `<FLUX_RESULTS>/figures/` |
| `GENE_ANNOT` | GENCODE v19 gene annotation TSV, for HGNC symbol labels | — (required by the scripts that label genes) |

```
GENE_ANNOT=/path/to/gencode.v19.genes.tsv \
FLUX_RESULTS=/path/to/results_findr_py \
python analysis/figure_scripts/plot_main_readouts.py
```

Figures are Arial, 300 DPI.

## Core analyses

| Script | What it does |
|--------|--------------|
| `core/build_flux_network.py` | builds one tissue's disease-regulatory flux network (node roles + flux edges) |
| `core/flux_cardiovascular_test.py` | cardiovascular functional-enrichment test and the curated gene universe |
| `core/fetch_omnipath_interactions.py` | downloads the OmniPath reference interaction network reproducibly |

## Main-text figures

| Script | Figure | What it shows |
|--------|--------|---------------|
| `figure_scripts/plot_prediction_summary.py` | `fig_alltargets_summary` | transcriptome-wide expression-prediction summary |
| `figure_scripts/plot_association_summary.py` | `fig_assoc_combined` | combined CAD transcriptome-wide association (Manhattan via `utils/`, needs `GENE_ANNOT`) |
| `figure_scripts/plot_flux_cascades.py` | `fig_flux_map_grn` | disease-regulatory flux map with the GRN cascade (needs `GENE_ANNOT`) |
| `figure_scripts/plot_flux_validation_findrpy.py` | `fig_flux_validation` | four-panel flux-map validation (OmniPath, sign-flip null, FinnGen replication) |
| `figure_scripts/plot_functional_axis.py` | `fig_functional_axis` | cardiovascular functional-relevance enrichment (needs `GENE_ANNOT`) |
| `figure_scripts/plot_main_readouts.py` | `fig_main_readouts` | six-panel disease-role readouts (needs `GENE_ANNOT`) |

## Supplementary figures

| Script | Figure | What it shows |
|--------|--------|---------------|
| `figure_scripts/plot_coexpression_density.py` | `fig_coexpression_density` | background co-expression versus reconstructed-network density per tissue |
| `figure_scripts/plot_prior_benchmark.py` | `fig_prior_benchmark` | predictive performance of the three shrinkage priors |
| `figure_scripts/plot_prior_pairplot.py` | `fig_prior_pairplot` | per-gene LOO-R² agreement between the three priors |
| `figure_scripts/plot_flux_lcc_network.py` | `fig_flux_lcc_<tissue>` (interactive HTML) | largest-connected-component flux map per tissue |

## Supplementary tables

The supplementary tables (S1–S14) are provided under `supplementary_tables/`, each with a
column-by-column data dictionary in `supplementary_tables/README.md`, and are archived on Zenodo.

## External validation

Both validators reuse the published predictor weights with no retraining and no access to the
individual-level STARNET data.

**Summary-statistic TWAS on a new trait** — `validation/twas_from_gwas.py` applies the weights to
any GWAS to obtain a per-gene, per-tissue TWAS statistic, decomposed into cis and trans
components. Align the GWAS to the ALT alleles first (`bin/harmonise_gwas.py`), then:

```
python analysis/validation/twas_from_gwas.py \
    --weights   <predictor weights CSV> \
    --sigma-g   <per-gene predicted-expression scale CSV> \
    --gwas      <trait>.aligned.tsv.gz \
    --out       <trait>.external_twas.csv
```

**Predictor portability in a new cohort** — `validation/validate_predictors_external.py` tests the
LOO-R² in an independent cohort with genotypes and measured expression from a matching tissue. It
predicts expression from dosage (cis and trans channels together) and reports the external squared
correlation against the measured expression, next to the reported `loo_r2`. No model is refit.

```
python analysis/validation/validate_predictors_external.py \
    --weights     <predictor weights CSV> \
    --loo         <per-gene LOO-R² CSV> \
    --dosage      <external dosage matrix.tsv.gz> \
    --expression  <external expression matrix.tsv.gz> \
    --tissue      AOR --out <cohort>.external_r2.csv --fig
```

## Shared helpers

Figure scripts import shared helpers from `utils/` (a repo-level package): `gene_labels.py`
(ENSG→HGNC symbols) and the Manhattan renderer (`plot_manhattan.py` + `manhattan_ggplot.R` +
`rscript_wrapper.sh`).
