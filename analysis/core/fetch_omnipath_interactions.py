#!/usr/bin/env python3
"""Download the OmniPath interaction network reproducibly, via the official
`omnipath` Python package (https://github.com/saezlab/omnipath).

It reproduces the six-dataset signalling + transcriptional network used throughout
this project -- literature-curated signalling (omnipath), unreferenced signalling
(pathwayextra, kinaseextra), ligand-receptor (ligrecextra) and transcription-factor
targets (collectri, dorothea levels A-C) -- and writes a single TSV plus a provenance
file so the exact build can be cited and re-obtained.

Each dataset is requested separately (the combined query overloads the server); the
results are concatenated. Downstream code deduplicates by gene-symbol pair, so any
interaction shared across datasets is counted once regardless.

The output schema matches what this project's validators read:
    source  target  source_genesymbol  target_genesymbol
    is_directed  is_stimulation  is_inhibition
    consensus_direction  consensus_stimulation  consensus_inhibition  sources  type

Usage:
    python analysis/core/fetch_omnipath_interactions.py \
        --out /cluster/projects/nn1015k/datasets/references/omnipath

Writes, under --out:
    omnipath_interactions.<YYYYMMDD>.tsv   the network
    omnipath_interactions.tsv              symlink/copy to the newest build
    README.txt                             provenance (package version, query,
                                           per-type row counts, sha256)
"""
import argparse
import hashlib
import os
import time
from datetime import date

import pandas as pd

# The six datasets, and the omnipath interaction class that serves each. dorothea is
# restricted to confidence levels A-C (highest three) to match the project's network.
DATASETS = ["omnipath", "pathwayextra", "kinaseextra", "ligrecextra", "collectri", "dorothea"]
DOROTHEA_LEVELS = ["A", "B", "C"]

# Columns kept, in fixed order.
COLUMNS = ["source", "target", "source_genesymbol", "target_genesymbol",
           "is_directed", "is_stimulation", "is_inhibition",
           "consensus_direction", "consensus_stimulation", "consensus_inhibition",
           "sources", "type"]


def fetch_one(op, name):
    """One dataset as a DataFrame with gene symbols and the `type` column."""
    common = dict(genesymbols=True, fields=["sources", "type"], license="academic")
    getters = {
        "omnipath":     lambda: op.interactions.OmniPath.get(**common),
        "pathwayextra": lambda: op.interactions.PathwayExtra.get(**common),
        "kinaseextra":  lambda: op.interactions.KinaseExtra.get(**common),
        "ligrecextra":  lambda: op.interactions.LigRecExtra.get(**common),
        "collectri":    lambda: op.interactions.CollecTRI.get(**common),
        "dorothea":     lambda: op.interactions.Dorothea.get(dorothea_levels=DOROTHEA_LEVELS, **common),
    }
    last = None
    for attempt in range(1, 7):
        try:
            return getters[name]()
        except Exception as e:                       # server is intermittently flaky
            last = e
            print(f"  {name}: attempt {attempt} failed ({type(e).__name__}); retrying", flush=True)
            time.sleep(20)
    raise RuntimeError(f"{name} failed after retries: {last}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="/cluster/projects/nn1015k/datasets/references/omnipath",
                    help="output directory")
    ap.add_argument("--stamp", default=None,
                    help="date stamp for the filename (default: today, YYYYMMDD)")
    args = ap.parse_args()
    os.makedirs(args.out, exist_ok=True)

    import omnipath as op
    op.options.progress_bar = False
    op.options.timeout = 180
    op.options.num_retries = 10

    frames, counts = [], {}
    for name in DATASETS:
        print(f"fetching {name} ...", flush=True)
        df = fetch_one(op, name)
        for c in COLUMNS:
            if c not in df.columns:
                df[c] = pd.NA
        frames.append(df[COLUMNS])
        counts[name] = len(df)
        print(f"  {name}: {len(df):,} interactions", flush=True)

    full = pd.concat(frames, ignore_index=True)
    stamp = args.stamp or date.today().strftime("%Y%m%d")
    out_tsv = os.path.join(args.out, f"omnipath_interactions.{stamp}.tsv")
    full.to_csv(out_tsv, sep="\t", index=False)

    # point the canonical filename at this build
    canonical = os.path.join(args.out, "omnipath_interactions.tsv")
    if os.path.islink(canonical) or os.path.exists(canonical):
        os.remove(canonical)
    os.symlink(os.path.basename(out_tsv), canonical)

    sha = hashlib.sha256(open(out_tsv, "rb").read()).hexdigest()
    by_type = full["type"].value_counts().to_dict()
    with open(os.path.join(args.out, "README.txt"), "w") as fh:
        fh.write(f"OmniPath interactions (omnipathdb.org), fetched {stamp} via the "
                 f"`omnipath` python package v{op.__version__}.\n")
        fh.write(f"Datasets: {', '.join(DATASETS)}  (dorothea levels {','.join(DOROTHEA_LEVELS)}); "
                 f"genesymbols=1; fields=sources,type; license=academic.\n")
        fh.write(f"Reproduce: python analysis/core/fetch_omnipath_interactions.py --out {args.out}\n\n")
        fh.write(f"File: {os.path.basename(out_tsv)}\n")
        fh.write(f"Rows: {len(full):,} total = " +
                 " + ".join(f"{v:,} {k}" for k, v in by_type.items()) + "\n")
        fh.write("Per-dataset rows (before cross-dataset union): " +
                 ", ".join(f"{k} {v:,}" for k, v in counts.items()) + "\n")
        fh.write(f"Columns: {', '.join(COLUMNS)}\n")
        fh.write(f"sha256({os.path.basename(out_tsv)}) = {sha}\n")

    print(f"\nwrote {out_tsv}\n  {len(full):,} interactions; sha256 {sha[:16]}...")
    print("  by type:", by_type)


if __name__ == "__main__":
    main()
