#!/bin/bash
# Run an Rscript with an R that provides ggplot2, ggrepel, dplyr and scales.
# Used by utils/plot_manhattan.py for the ggrepel Manhattan renderer.
#
# On an Lmod-based HPC this loads an R + Bioconductor module; elsewhere it falls back
# to whatever `Rscript` is already on PATH. Adapt the module line to your own R install.
if [ -f /cluster/installations/lmod/lmod/init/bash ]; then
    source /cluster/installations/lmod/lmod/init/bash
    module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2 2>/dev/null || true
fi
exec Rscript "$@"
