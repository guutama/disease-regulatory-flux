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

## For development

| software | version | role |
|----------|---------|------|
| pytest | 8.4.2 | runs the unit tests in `tests/` |

## Per-step requirements

| step | requirements |
|------|--------------|
| 1. harmonise inputs (`bin/harmonise_inputs.py`) | Python ≥ 3.9, standard library only |
| 2. select instruments (`bin/select_instruments.py`) | Python ≥ 3.9, standard library only |
