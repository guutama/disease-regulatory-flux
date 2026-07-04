# Analysis code

Scripts that produce the figures and tables in the disease-regulatory-flux manuscript, kept
separate from the reproducible pipeline (`bin/`, `modules/`, `flux.nf`). Each script is the code
as run, with a module docstring saying what it makes, which inputs it reads, and how to run it.

```
analysis/
  figures/   one script per manuscript figure (main text and supplementary)
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

## Shared helpers

Figure scripts import shared helpers from `utils/` (a repo-level package): `gene_labels.py`
(ENSG→HGNC symbols), `figure_data.py` (shared data layer and style), and the Manhattan renderer
(`plot_manhattan.py` + `manhattan_ggplot.R` + `rscript_wrapper.sh`).
