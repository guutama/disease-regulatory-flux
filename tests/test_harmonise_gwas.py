"""Unit tests for bin/harmonise_gwas.py (Stage 8: GWAS harmonisation).

The pipeline genotype counts the ALT allele, so every GWAS effect must be re-expressed
as the effect of the pipeline's ALT allele before the TWAS. We test:
  * complement / strand-ambiguity helpers,
  * the allele alignment rule (direct match, strand flip, ambiguity drop, mismatch drop),
  * the end-to-end run: matching a GWAS to the variant universe by chrom/pos, aligning the
    effect (and z) to ALT, dropping ambiguous and unmatched variants, value pass-through.
"""
import gzip
import importlib.util
from pathlib import Path

import pytest

_SPEC = importlib.util.spec_from_file_location(
    "harmonise_gwas", Path(__file__).resolve().parents[1] / "bin" / "harmonise_gwas.py")
hg = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(hg)


def _write(path: Path, rows: list[list[str]]) -> None:
    text = "\n".join("\t".join(map(str, r)) for r in rows) + "\n"
    if str(path).endswith(".gz"):
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)


# -------------------------------------------------------------------- helper-function tests
def test_complement():
    assert hg.complement("A") == "T"
    assert hg.complement("t") == "A"
    assert hg.complement("C") == "G"
    assert hg.complement("G") == "C"
    assert hg.complement("N") is None          # not a simple base


def test_strand_ambiguous():
    assert hg.is_ambiguous("A", "T")
    assert hg.is_ambiguous("C", "G")
    assert hg.is_ambiguous("g", "c")
    assert not hg.is_ambiguous("A", "G")
    assert not hg.is_ambiguous("A", "C")


def test_align_sign_direct():
    # effect allele is ALT -> +1 ; effect allele is REF -> -1
    assert hg.align_sign("G", "A", "A", "G") == 1
    assert hg.align_sign("A", "G", "A", "G") == -1


def test_align_sign_strand_flip():
    # GWAS reports the complementary strand (C/T) of a REF/ALT = G/A SNP
    assert hg.align_sign("T", "C", "G", "A") == 1     # T complements A (=ALT) -> +1
    assert hg.align_sign("C", "T", "G", "A") == -1    # C complements G (=REF) -> -1


def test_align_sign_ambiguous_and_mismatch():
    assert hg.align_sign("A", "T", "A", "T") is None  # strand-ambiguous: cannot resolve
    assert hg.align_sign("A", "C", "A", "G") is None  # alleles do not match REF/ALT


# ------------------------------------------------------------------------- end-to-end tests
# Variant universe (e.g. a genotype_012 matrix): only the id + allele columns are read.
UNIVERSE = [["variant_id", "chromosome", "position", "ref", "alt", "S1", "S2"],
            ["v_alt",  "1", "100", "A", "G", 0, 1],   # GWAS effect on ALT (G)
            ["v_ref",  "1", "200", "C", "T", 1, 2],   # GWAS effect on REF (C) -> flip
            ["v_flip", "2", "300", "G", "A", 0, 0],   # GWAS on complementary strand
            ["v_amb",  "2", "400", "A", "T", 1, 1],   # strand-ambiguous -> dropped
            ["v_miss", "3", "500", "C", "G", 2, 0]]   # no GWAS row -> dropped

# GWAS: effect_allele / other_allele, beta, se, p, n.
GWAS = [["chromosome", "position", "effect_allele", "other_allele", "beta", "se", "pvalue", "n"],
        ["1", "100", "G", "A", "0.20", "0.10", "1e-5", "1000"],   # v_alt: ea=ALT -> z=+2
        ["1", "200", "C", "T", "0.30", "0.10", "1e-8", "1000"],   # v_ref: ea=REF -> z=-3
        ["2", "300", "T", "C", "0.40", "0.10", "2e-4", "900"],    # v_flip: T~ALT -> z=+4
        ["2", "400", "A", "T", "0.50", "0.10", "3e-3", "900"]]    # v_amb: ambiguous -> drop


def run(d: Path, gwas=GWAS, universe=UNIVERSE, extra=None):
    _write(d / "gwas.tsv", gwas)
    _write(d / "uni.tsv.gz", universe)
    out = d / "CADtest.gwas.tsv.gz"
    qc = d / "CADtest.qc.tsv"
    argv = ["--gwas", str(d / "gwas.tsv"), "--variants", str(d / "uni.tsv.gz"),
            "--trait", "CADtest", "--out", str(out), "--qc", str(qc)]
    if extra:
        argv += extra
    hg.main(argv)
    rows = [ln.split("\t") for ln in gzip.decompress(out.read_bytes()).decode().splitlines()]
    qc_d = {k: v for k, v in (ln.split("\t") for ln in qc.read_text().splitlines()[1:])}
    return rows, qc_d


def test_output_header_and_kept_variants(tmp_path):
    rows, qc = run(tmp_path)
    assert rows[0] == ["variant_id", "z", "beta", "se", "p_value", "n"]
    kept = {r[0] for r in rows[1:]}
    # v_alt, v_ref, v_flip kept; v_amb (ambiguous) and v_miss (no GWAS) dropped
    assert kept == {"v_alt", "v_ref", "v_flip"}


def test_alt_alignment_signs(tmp_path):
    rows, _ = run(tmp_path)
    z = {r[0]: float(r[1]) for r in rows[1:]}
    assert z["v_alt"] == pytest.approx(2.0)    # ea=ALT: +beta/se
    assert z["v_ref"] == pytest.approx(-3.0)   # ea=REF: -beta/se
    assert z["v_flip"] == pytest.approx(4.0)   # strand flip, effect on ALT: +beta/se


def test_value_passthrough(tmp_path):
    rows, _ = run(tmp_path)
    r = {x[0]: x for x in rows[1:]}
    # v_ref: beta sign-flipped to ALT (-0.30), se/p/n unchanged
    assert float(r["v_ref"][2]) == pytest.approx(-0.30)
    assert r["v_ref"][3:] == ["0.10", "1e-8", "1000"]


def test_qc_counts(tmp_path):
    _, qc = run(tmp_path)
    assert qc["trait"] == "CADtest"
    assert qc["variants_in_universe"] == "5"
    assert qc["aligned"] == "3"
    assert qc["dropped_ambiguous"] == "1"
    assert qc["dropped_unmatched"] == "1"
