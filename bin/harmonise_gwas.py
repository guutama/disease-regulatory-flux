#!/usr/bin/env python3
"""Stage 8 -- harmonise a GWAS summary-statistics file to the pipeline's ALT allele.

The pipeline genotype counts the ALT allele and the predictor weights are in ALT-dosage
units, so every GWAS effect must be re-expressed as the effect of the pipeline's ALT
allele before the TWAS association. For each variant in the pipeline's variant universe
(the harmonised genotype / genotype_012 matrices, which carry chromosome, position, ref
and alt), this script finds the matching GWAS record by chromosome and position, checks
the alleles, and writes the GWAS z-score aligned to ALT:

  * effect on ALT        -> z = +beta/se
  * effect on REF        -> z = -beta/se   (sign flipped)
  * complementary strand -> resolved by complementing the GWAS alleles
  * strand-ambiguous SNPs (A/T, C/G) -> dropped (the strand cannot be resolved)
  * alleles that do not match REF/ALT, or have no GWAS record -> dropped

Inputs
  --gwas       raw GWAS summary statistics (.tsv/.tsv.gz), one row per variant
  --variants   one or more files carrying variant_id, chromosome, position, ref, alt
               (e.g. the per-tissue genotype_012 matrices); only those columns are read
  --trait      trait label (recorded in the QC report)
  --out        output path (.tsv.gz)
  --qc         QC report path

Output (--out): one row per aligned variant, with z and beta on the ALT scale
  variant_id  z  beta  se  p_value  n
"""
from __future__ import annotations

import argparse
import gzip
from pathlib import Path

_BASES = {"A", "C", "G", "T"}
_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _open(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _split(line: str, delim: str) -> list[str]:
    return line.rstrip("\n").split(delim)


# ----------------------------------------------------------------------------- allele logic
def complement(allele: str) -> str | None:
    """Watson-Crick complement of a single base (case-insensitive), or None if not A/C/G/T."""
    return _COMP.get(allele.upper())


def is_ambiguous(a: str, b: str) -> bool:
    """True if the two alleles are strand-ambiguous (A/T or C/G), so strand is unresolvable."""
    return {a.upper(), b.upper()} in ({"A", "T"}, {"C", "G"})


def alleles_match(ea: str, oa: str, ref: str, alt: str) -> bool:
    """True if the GWAS alleles {ea, oa} are REF/ALT, directly or on the complementary strand."""
    ea, oa, ref, alt = ea.upper(), oa.upper(), ref.upper(), alt.upper()
    if not {ea, oa, ref, alt} <= _BASES:
        return False
    if {ea, oa} == {ref, alt}:
        return True
    ce, co = complement(ea), complement(oa)
    return ce is not None and co is not None and {ce, co} == {ref, alt}


def align_sign(ea: str, oa: str, ref: str, alt: str) -> int | None:
    """Sign to put a GWAS effect on the ALT allele: +1 (effect on ALT), -1 (effect on REF),
    or None when the SNP is strand-ambiguous or the alleles do not match REF/ALT."""
    if not alleles_match(ea, oa, ref, alt):
        return None
    if is_ambiguous(ref, alt):
        return None
    ea, ref, alt = ea.upper(), ref.upper(), alt.upper()
    if ea == alt:
        return 1            # effect already on ALT
    if ea == ref:
        return -1           # effect on REF -> flip
    # complementary strand: the effect allele complements ALT (+1) or REF (-1)
    return 1 if complement(ea) == alt else -1


def norm_chrom(c: str) -> str:
    """Drop a leading 'chr' so '1' and 'chr1' match."""
    c = c.strip()
    return c[3:] if c.lower().startswith("chr") else c


def norm_pos(p: str) -> str:
    """Normalise a base-pair position (drop leading zeros) so '100' and '0100' match."""
    p = p.strip()
    return str(int(p)) if p.isdigit() else p


# -------------------------------------------------------------------------------- universe
def load_universe(paths: list[str], delim: str, variant_col: str):
    """Read variant_id/chromosome/position/ref/alt from each file. Returns (positions, n)
    where positions[(chrom, pos)] is a list of (variant_id, ref, alt). The first-seen
    variant id wins on duplicates; n is the number of distinct variants."""
    positions: dict[tuple[str, str], list[tuple[str, str, str]]] = {}
    seen: set[str] = set()
    n = 0
    for path in paths:
        with _open(path) as fh:
            header = _split(fh.readline(), delim)
            idx = {c: i for i, c in enumerate(header)}
            for col in (variant_col, "chromosome", "position", "ref", "alt"):
                if col not in idx:
                    raise SystemExit(
                        f"The variants file {path} must carry a '{col}' column; "
                        f"it has {header[:6]}.")
            vi, ci, pi, ri, ai = (idx[variant_col], idx["chromosome"],
                                  idx["position"], idx["ref"], idx["alt"])
            for line in fh:
                f = _split(line, delim)
                if not f or f == [""]:
                    continue
                vid = f[vi]
                if vid in seen:
                    continue
                seen.add(vid)
                n += 1
                key = (norm_chrom(f[ci]), norm_pos(f[pi]))
                positions.setdefault(key, []).append((vid, f[ri], f[ai]))
    return positions, n


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gwas", required=True, help="input: raw GWAS summary statistics")
    p.add_argument("--variants", required=True, nargs="+",
                   help="input: one or more files with variant_id/chromosome/position/ref/alt "
                        "(e.g. the genotype_012 matrices)")
    p.add_argument("--trait", required=True, help="trait label (recorded in the QC report)")
    p.add_argument("--out", required=True, help="output: ALT-aligned GWAS (.tsv.gz)")
    p.add_argument("--qc", required=True, help="output: QC report")
    p.add_argument("--delimiter", default="\t", help="field delimiter (default: tab)")
    p.add_argument("--variant-col", default="variant_id",
                   help="variant-id column in the variants files (default: variant_id)")
    p.add_argument("--chrom-col", default="chromosome", help="GWAS chromosome column")
    p.add_argument("--pos-col", default="position", help="GWAS position column")
    p.add_argument("--effect-allele-col", default="effect_allele",
                   help="GWAS effect-allele column (the allele the effect refers to)")
    p.add_argument("--other-allele-col", default="other_allele", help="GWAS other-allele column")
    p.add_argument("--beta-col", default="beta", help="GWAS effect-size column")
    p.add_argument("--se-col", default="se", help="GWAS standard-error column")
    p.add_argument("--p-col", default="pvalue", help="GWAS p-value column")
    p.add_argument("--n-col", default="n", help="GWAS sample-size column (optional)")
    a = p.parse_args(argv)
    if a.delimiter == "\\t":
        a.delimiter = "\t"
    return a


def main(argv=None) -> None:
    args = parse_args(argv)
    d = args.delimiter

    positions, n_universe = load_universe(args.variants, d, args.variant_col)

    # ---- stream the GWAS, keeping only records at a universe position ------------------
    gwas_at: dict[tuple[str, str], list[tuple]] = {}
    with _open(args.gwas) as fh:
        header = _split(fh.readline(), d)
        idx = {c: i for i, c in enumerate(header)}
        for col in (args.chrom_col, args.pos_col, args.effect_allele_col,
                    args.other_allele_col, args.beta_col, args.se_col, args.p_col):
            if col not in idx:
                raise SystemExit(f"The GWAS file must carry a '{col}' column; it has {header}.")
        ni = idx.get(args.n_col)
        ci, pi = idx[args.chrom_col], idx[args.pos_col]
        eai, oai = idx[args.effect_allele_col], idx[args.other_allele_col]
        bi, si, ppi = idx[args.beta_col], idx[args.se_col], idx[args.p_col]
        for line in fh:
            f = _split(line, d)
            if not f or f == [""]:
                continue
            key = (norm_chrom(f[ci]), norm_pos(f[pi]))
            if key not in positions:
                continue
            n_val = f[ni] if ni is not None and ni < len(f) else ""
            gwas_at.setdefault(key, []).append(
                (f[eai], f[oai], f[bi], f[si], f[ppi], n_val))

    # ---- align each universe variant to its GWAS record -------------------------------
    n_aligned = n_ambiguous = n_unmatched = 0
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with _open(out, "wt") as fout:
        fout.write(d.join(["variant_id", "z", "beta", "se", "p_value", "n"]) + "\n")
        for key, variants in positions.items():
            records = gwas_at.get(key, [])
            for vid, ref, alt in variants:
                rec = next((r for r in records if alleles_match(r[0], r[1], ref, alt)), None)
                if rec is None:
                    n_unmatched += 1
                    continue
                if is_ambiguous(ref, alt):
                    n_ambiguous += 1
                    continue
                ea, oa, beta, se, pval, n_val = rec
                sign = align_sign(ea, oa, ref, alt)
                z = sign * float(beta) / float(se)
                beta_alt = sign * float(beta)
                fout.write(d.join([vid, f"{z:.6g}", f"{beta_alt:.6g}",
                                   se, pval, n_val]) + "\n")
                n_aligned += 1

    with open(args.qc, "wt") as fh:
        fh.write("metric\tvalue\n")
        for k, v in (("trait", args.trait),
                     ("variants_in_universe", n_universe),
                     ("aligned", n_aligned),
                     ("dropped_ambiguous", n_ambiguous),
                     ("dropped_unmatched", n_unmatched)):
            fh.write(f"{k}\t{v}\n")
    print(f"[harmonise_gwas] trait={args.trait} universe={n_universe} aligned={n_aligned} "
          f"ambiguous={n_ambiguous} unmatched={n_unmatched}")


if __name__ == "__main__":
    main()
