# Requirements

Software needed to run the pipeline. The listed versions are what it is developed and
tested against; the same or newer should work. How you install them — directly, with conda,
or in a container — is up to you.

## To run the pipeline

| software | version | role |
|----------|---------|------|
| Nextflow | ≥ 24.04.2 (tested with 24.04.2) | workflow engine |
| Java     | 11 or newer | required by Nextflow |
| Python   | ≥ 3.9 (tested with 3.9.21) | runs the pipeline scripts |
| NumPy    | ≥ 1.21 | LD-filtering and expression-model steps |
| NumPyro + JAX | NumPyro ≥ 0.13 | MCMC fitting of the expression models |
| ArviZ    | 0.17 (with SciPy < 1.13), or a build matching your SciPy | PSIS leave-one-out scoring of the expression models |
| Julia    | 1.11.3 | runs the GRN reconstruction step (BioFindr) |
| BioFindr (+ CSV, DataFrames, Graphs) | pinned in `julia_project/` | causal network inference |

ArviZ 0.17 imports `scipy.signal.gaussian`, which moved to `scipy.signal.windows` in SciPy
1.13; install ArviZ against SciPy < 1.13 (or a newer ArviZ matched to your SciPy) so it loads.

The Julia packages are pinned in `julia_project/` (`Project.toml` / `Manifest.toml`) and
installed once with `julia_project/install_biofindr.jl`. Julia 1.11.3 is required: earlier
1.10.x releases fail to load the pinned BioFindr dependency tree.

## For development

| software | version | role |
|----------|---------|------|
| pytest | 8.4.2 | runs the unit tests in `tests/` |

## Per-step requirements

| step | requirements |
|------|--------------|
| 1. harmonise inputs (`bin/harmonise_inputs.py`) | Python ≥ 3.9, standard library only |
| 2. select instruments (`bin/select_instruments.py`) | Python ≥ 3.9, standard library only |
| 3. reconstruct GRN (`bin/reconstruct_grn.jl`) | Julia 1.11.3 + BioFindr (pinned in `julia_project/`) |
| 4. select cis features (`bin/select_cis_features.py`) | Python ≥ 3.9, standard library only |
| 5. LD-prune (`bin/ld_prune.py`) | Python ≥ 3.9 + NumPy |
| 6. trans features (`bin/trans_features.py`) | Python ≥ 3.9, standard library only |
| 7. expression models (`bin/fit_expression_models.py`) | Python ≥ 3.9 + NumPy, NumPyro + JAX (MCMC), ArviZ (PSIS-LOO) |
| 8. harmonise GWAS (`bin/harmonise_gwas.py`) | Python ≥ 3.9, standard library only |
| 9. TWAS association (`bin/twas_association.py`) | Python ≥ 3.9 + NumPy |
| 10. flux map (`bin/build_flux_map.py`) | Python ≥ 3.9 + NumPy |
