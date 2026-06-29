# Requirements

Software needed to run the pipeline. The listed versions are what it is developed and
tested against; the same or newer should work. How you install them — directly, with conda,
or in a container — is up to you.

## To run the pipeline

| software | version | role |
|----------|---------|------|
| Nextflow | ≥ 24.04.2 (tested with 24.04.2) | workflow engine |
| Java     | 11 or newer | required by Nextflow |
| Python   | ≥ 3.9 (tested with 3.9.21) | runs the pipeline scripts; standard library only |
| Julia    | 1.11.3 | runs the GRN reconstruction step (BioFindr) |
| BioFindr (+ CSV, DataFrames, Graphs) | pinned in `julia_project/` | causal network inference |

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
