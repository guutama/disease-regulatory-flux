#!/usr/bin/env python3
"""
Per-tissue independent random train/test split of a tissue's samples.

The samples to split can be given two ways:

  --genotype    one or more 0/1/2 genotype matrices whose header is
                variant_id, chromosome, position, ref, alt, then one column per
                sample; the sample ids are those columns. This is the flux.nf
                path (the samples never need to be listed separately).
  --samples-tsv one or more <tissue>.samples.tsv audit maps (columns
                genotype_id + kept); the kept=1 ids are split.

Writes per tissue:
    <tissue>.train.samples    one sample id per line
    <tissue>.test.samples     one sample id per line (empty when the whole
                              cohort is used for training; see below)
plus
    <tissue>.split_qc.json    (single-tissue run) or split_qc.json (multi):
                              per-tissue counts + manifest (seed, fraction).

Train fraction:
  --train-fraction is in (0, 1]. A fraction below 1 assigns that share of each
  tissue's samples to train and the rest to test, so downstream steps evaluate
  predictors on the held-out test fold. A fraction of exactly 1 assigns every
  sample to train and writes an EMPTY test file: the whole cohort is used for
  both fitting and (leave-one-out) scoring, which is the intended behaviour when
  the cohort is too small to spare a held-out fold.

Reproducibility:
  - One numpy RandomState per tissue, derived from (split_seed, tissue) via a
    deterministic hash, so:
      - Re-running the same seed yields identical splits, AND
      - Changing the seed changes ALL tissues' splits (no accidental coupling).
  - Tissues are processed in alphabetical order, regardless of input order.

The split is OVER THE TISSUE'S SAMPLES ONLY. A sample may end up in train for
one tissue and test for another -- that is the "per-tissue independent"
semantics. If cross-tissue consistency is needed later, switch to a
shared-split variant.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import signal
import sys
from pathlib import Path

import numpy as np

try:
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except (AttributeError, ValueError):
    pass

# Genotype-matrix filename suffixes, longest first, so <tissue> is recovered cleanly.
GENO_SUFFIXES = (".genotype_012.tsv.gz", ".genotype_012.tsv", ".genotype.tsv.gz",
                 ".genotype.tsv", ".tsv.gz", ".tsv")
# Fixed leading columns of a genotype matrix before the per-sample columns.
GENO_LEADING_COLS = 5


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--genotype", nargs="+",
                     help="one or more genotype matrices; sample ids are the columns after "
                          "variant_id, chromosome, position, ref, alt")
    src.add_argument("--samples-tsv", nargs="+",
                     help="one or more <tissue>.samples.tsv files (columns genotype_id + kept)")
    p.add_argument("--tissue", default=None,
                   help="tissue label for a single --genotype file (default: from the filename)")
    p.add_argument("--train-fraction", type=float, required=True,
                   help="fraction (0,1] of each tissue's samples assigned to train; "
                        "1 uses the whole cohort for train and writes an empty test file")
    p.add_argument("--split-seed", type=int, required=True,
                   help="base seed; combined with tissue name -> per-tissue RNG")
    p.add_argument("--outdir", default=".")
    return p.parse_args()


def _open(path: str):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_kept_samples(samples_tsv: str) -> list[str]:
    """genotype_id of every kept=1 row, in file order."""
    out: list[str] = []
    with open(samples_tsv) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_gid  = header.index("genotype_id")
            i_kept = header.index("kept")
        except ValueError as e:
            sys.exit(f"ERROR: {samples_tsv} missing column: {e}")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[i_kept] == "1":
                out.append(f[i_gid])
    return out


def load_genotype_samples(genotype: str) -> list[str]:
    """Sample ids = the genotype matrix's header columns after the fixed leading columns."""
    with _open(genotype) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    if len(header) <= GENO_LEADING_COLS:
        sys.exit(f"ERROR: {genotype} has no sample columns after the first {GENO_LEADING_COLS}")
    return header[GENO_LEADING_COLS:]


def tissue_from_geno(name: str) -> str:
    for suf in GENO_SUFFIXES:
        if name.endswith(suf):
            return name[:-len(suf)]
    return name.split(".")[0]


def collect_tissues(args) -> list[tuple[str, list[str]]]:
    """Return (tissue, samples) for every input, from whichever source was given."""
    pairs: list[tuple[str, list[str]]] = []
    if args.samples_tsv:
        for f in args.samples_tsv:
            name = Path(f).name
            if not name.endswith(".samples.tsv"):
                sys.exit(f"ERROR: input filename must end with .samples.tsv: {f}")
            pairs.append((name[:-len(".samples.tsv")], load_kept_samples(f)))
    else:
        if args.tissue and len(args.genotype) > 1:
            sys.exit("ERROR: --tissue is only valid with a single --genotype file")
        for f in args.genotype:
            tissue = args.tissue if args.tissue else tissue_from_geno(Path(f).name)
            pairs.append((tissue, load_genotype_samples(f)))
    pairs.sort(key=lambda p: p[0])   # alphabetical for reproducibility
    return pairs


def tissue_seed(base_seed: int, tissue: str) -> int:
    """Deterministic per-tissue seed in [0, 2^31). Stable across Python runs
    (uses hashlib, not the salted built-in hash())."""
    h = hashlib.sha256(f"{base_seed}|{tissue}".encode()).digest()
    return int.from_bytes(h[:4], "big") & 0x7FFFFFFF


def main() -> None:
    args = parse_args()
    if not (0.0 < args.train_fraction <= 1.0):
        sys.exit(f"ERROR: --train-fraction must be in (0,1]; got {args.train_fraction}")
    full_cohort = args.train_fraction == 1.0
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    tissues = collect_tissues(args)
    qc_tissues = []
    for tissue, samples in tissues:
        min_needed = 1 if full_cohort else 2
        if len(samples) < min_needed:
            sys.exit(f"ERROR: [{tissue}] only {len(samples)} samples; too few to split")

        seed = tissue_seed(args.split_seed, tissue)
        if full_cohort:
            # Whole cohort to train; empty test file signals leave-one-out scoring downstream.
            train_ids = list(samples)
            test_ids: list[str] = []
        else:
            n_train = int(round(args.train_fraction * len(samples)))
            # Enforce at least 1 in each split so degenerate edge cases don't silently
            # produce an empty file.
            n_train = max(1, min(len(samples) - 1, n_train))
            rng = np.random.default_rng(seed)
            order = np.arange(len(samples))
            rng.shuffle(order)
            train_idx = sorted(order[:n_train].tolist())
            test_idx  = sorted(order[n_train:].tolist())
            train_ids = [samples[i] for i in train_idx]
            test_ids  = [samples[i] for i in test_idx]

        (outdir / f"{tissue}.train.samples").write_text("\n".join(train_ids) + "\n")
        # An empty test split writes an empty file (no stray blank line).
        (outdir / f"{tissue}.test.samples").write_text(
            ("\n".join(test_ids) + "\n") if test_ids else "")

        qc_tissues.append({
            "tissue": tissue,
            "n_total": len(samples),
            "n_train": len(train_ids),
            "n_test":  len(test_ids),
            "tissue_seed": seed,
        })
        sys.stderr.write(
            f"[{tissue}] n={len(samples)} -> train={len(train_ids)} "
            f"test={len(test_ids)} (seed={seed})\n"
        )

    manifest = {
        "base_seed":      args.split_seed,
        "train_fraction": args.train_fraction,
        "n_tissues":      len(qc_tissues),
        "tissues":        qc_tissues,
    }
    # Single-tissue runs (the flux.nf per-tissue task) name the manifest per tissue so
    # published files do not collide; a multi-tissue run keeps the shared name.
    qc_name = f"{qc_tissues[0]['tissue']}.split_qc.json" if len(qc_tissues) == 1 else "split_qc.json"
    (outdir / qc_name).write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
