# Analysis code

Scripts that produce the figures and tables in the disease-regulatory-flux manuscript, kept
separate from the pipeline (`bin/`, `modules/`). Each script is the code as run, with a docstring
saying what it makes and which inputs it reads.

## Data availability

STARNET is individual-level, controlled-access data (dbGaP/EGA); raw data requires data-access
approval and is not included in this repository.

## Paths

Result locations resolve under the `FLUX_RESULTS` environment variable, defaulting to the repo's
own `results_cv/`. Override it to point elsewhere:

```
FLUX_RESULTS=/path/to/results python analysis/figures/plot_alltargets_summary.py
```

## Figures and tables

| Script | Manuscript item |
|--------|-----------------|
| `figures/plot_alltargets_summary.py` | Fig. `fig:pred` — transcriptome-wide prediction summary |
