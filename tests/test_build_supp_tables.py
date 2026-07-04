"""Unit tests for the Supplementary Table S9 loaders in bin/build_supp_tables.py.

S9 lists, for every predictable gene, the per-SNP weights of the gene's *selected*
predictor configuration, annotated with dbSNP coordinates and alleles taken from a
variant map. These tests cover the two functions that build it: load_snp_coords (the
streamed variant-map lookup) and load_weights (the per-tissue selection + annotation).
"""
import importlib.util
from pathlib import Path

import pandas as pd
import pytest

_SPEC = importlib.util.spec_from_file_location(
    "build_supp_tables",
    Path(__file__).resolve().parents[1] / "analysis" / "tables" / "build_supp_tables.py")
bst = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(bst)


def _write_tsv(path, rows, gzip_it=False):
    df = pd.DataFrame(rows)
    df.to_csv(path, sep="\t", index=False, compression="gzip" if gzip_it else None)


@pytest.mark.parametrize("gz", [False, True])
def test_load_snp_coords_keeps_only_requested(tmp_path, gz):
    p = tmp_path / ("map.tsv.gz" if gz else "map.tsv")
    _write_tsv(p, [
        {"original_id": "rsA", "chromosome": "1", "position": "100", "ref": "A", "alt": "G"},
        {"original_id": "rsB", "chromosome": "2", "position": "200", "ref": "C", "alt": "T"},
        {"original_id": "rsC", "chromosome": "3", "position": "300", "ref": "G", "alt": "A"},
    ], gzip_it=gz)
    coords = bst.load_snp_coords(str(p), "original_id", {"rsA", "rsC"})
    assert set(coords) == {"rsA", "rsC"}          # rsB (not requested) dropped
    assert coords["rsA"] == ("1", "100", "A", "G")
    assert coords["rsC"] == ("3", "300", "G", "A")


def _tissue_fixtures(d, monkeypatch):
    """One tissue T1 with three genes exercising every selection branch."""
    monkeypatch.setattr(bst, "TISSUES", ["T1"])
    monkeypatch.setattr(bst, "MODELS", str(d))
    monkeypatch.setattr(bst, "WEIGHTS", str(d / "{t}.horseshoe_weights.tsv.gz"))
    monkeypatch.setattr(bst, "FEATURES", str(d / "{t}.trans_features.tsv.gz"))
    monkeypatch.setattr(bst, "SNP_SD", str(d / "{t}.snp_sd.tsv.gz"))

    _write_tsv(d / "T1.horseshoe_metrics.tsv.gz", [
        {"gene_id": "g1", "best_config": "cis_only", "predictable": "True"},
        {"gene_id": "g2", "best_config": "trans_only", "predictable": "False"},
        {"gene_id": "g3", "best_config": "cis_trans", "predictable": "True"},
    ], gzip_it=True)

    _write_tsv(d / "T1.horseshoe_weights.tsv.gz", [
        # g1 selected cis_only -> only rs1/rs2 survive, its trans_only rs3 is dropped
        {"gene_id": "g1", "method": "horseshoe", "config": "cis_only",  "channel": "cis",   "rs_id": "rs1", "weight": 0.5},
        {"gene_id": "g1", "method": "horseshoe", "config": "cis_only",  "channel": "cis",   "rs_id": "rs2", "weight": -0.2},
        {"gene_id": "g1", "method": "horseshoe", "config": "trans_only", "channel": "trans", "rs_id": "rs3", "weight": 0.9},
        # g2 is not predictable -> excluded entirely
        {"gene_id": "g2", "method": "horseshoe", "config": "trans_only", "channel": "trans", "rs_id": "rs4", "weight": 0.1},
        # g3 selected cis_trans -> cis rs5 (mapped) + two trans SNPs survive; its cis_only rs6 dropped.
        # rs_absent is a trans SNP tagging one regulator but missing from the variant map;
        # rs7 is a trans SNP tagging two regulators (';'-joined).
        {"gene_id": "g3", "method": "horseshoe", "config": "cis_trans", "channel": "cis",   "rs_id": "rs5",      "weight": 0.3},
        {"gene_id": "g3", "method": "horseshoe", "config": "cis_trans", "channel": "trans", "rs_id": "rs_absent", "weight": 0.4},
        {"gene_id": "g3", "method": "horseshoe", "config": "cis_trans", "channel": "trans", "rs_id": "rs7",      "weight": 0.6},
        {"gene_id": "g3", "method": "horseshoe", "config": "cis_only",  "channel": "cis",   "rs_id": "rs6", "weight": 0.7},
    ], gzip_it=True)

    _write_tsv(d / "T1.trans_features.tsv.gz", [
        {"gene_id": "g3", "hop": 1, "source_gene_id": "gR",  "rs_id": "rs_absent"},
        {"gene_id": "g3", "hop": 2, "source_gene_id": "gR2", "rs_id": "rs7"},
        {"gene_id": "g3", "hop": 1, "source_gene_id": "gR",  "rs_id": "rs7"},
    ], gzip_it=True)

    # snp_sd omits rs_absent, so that SNP has both a missing allele and a blank sd.
    _write_tsv(d / "T1.snp_sd.tsv.gz", [
        {"rs_id": "rs1", "sd": 1.0},
        {"rs_id": "rs2", "sd": 0.5},
        {"rs_id": "rs5", "sd": 0.8},
        {"rs_id": "rs7", "sd": 0.7},
    ], gzip_it=True)

    vm = d / "map.tsv.gz"
    _write_tsv(vm, [
        {"original_id": "rs1", "chromosome": "1", "position": "10", "ref": "A", "alt": "G"},
        {"original_id": "rs2", "chromosome": "1", "position": "20", "ref": "C", "alt": "T"},
        {"original_id": "rs5", "chromosome": "5", "position": "50", "ref": "G", "alt": "A"},
        {"original_id": "rs7", "chromosome": "7", "position": "70", "ref": "T", "alt": "A"},
    ], gzip_it=True)
    return vm


def test_load_weights_selection_and_annotation(tmp_path, monkeypatch):
    vm = _tissue_fixtures(tmp_path, monkeypatch)
    mapping = {"g1": "GENE1", "g3": "GENE3", "gR": "REG", "gR2": "REG2"}

    s9 = bst.load_weights(mapping, str(vm), "original_id")

    assert list(s9.columns) == bst.S9_COLS
    # g2 (not predictable) never appears
    assert set(s9["gene_id"]) == {"g1", "g3"}
    # exactly the selected configuration per gene, nothing else
    assert set(s9.loc[s9.gene_id == "g1", "rs_id"]) == {"rs1", "rs2"}
    assert set(s9.loc[s9.gene_id == "g3", "rs_id"]) == {"rs5", "rs_absent", "rs7"}
    assert s9.groupby("gene_id")["config"].nunique().max() == 1
    # alleles annotated from the map
    rs1 = s9[s9.rs_id == "rs1"].iloc[0]
    assert (rs1["chromosome"], rs1["position"], rs1["ref"], rs1["alt"]) == ("1", "10", "A", "G")
    # a SNP absent from the map is retained with blank alleles and blank sd
    absent = s9[s9.rs_id == "rs_absent"].iloc[0]
    assert pd.isna(absent["chromosome"]) and pd.isna(absent["alt"])
    assert pd.isna(absent["sd"])
    # per-SNP sd is joined from the snp_sd table
    assert s9[s9.rs_id == "rs1"].iloc[0]["sd"] == 1.0
    # symbols resolved from the mapping
    assert set(s9.loc[s9.gene_id == "g1", "gene"]) == {"GENE1"}

    # cis SNPs carry no regulator; the channel marks them
    rs5 = s9[s9.rs_id == "rs5"].iloc[0]
    assert rs5["channel"] == "cis" and rs5["source_gene"] == "" and pd.isna(rs5["hop"])
    # a trans SNP names its regulator and hop
    absent_reg = s9[s9.rs_id == "rs_absent"].iloc[0]
    assert absent_reg["channel"] == "trans"
    assert absent_reg["source_gene"] == "REG" and absent_reg["source_gene_id"] == "gR"
    assert absent_reg["hop"] == 1
    # a trans SNP tagging two regulators lists both (sorted, ';'-joined) with the minimum hop
    rs7 = s9[s9.rs_id == "rs7"].iloc[0]
    assert rs7["source_gene_id"] == "gR;gR2" and rs7["source_gene"] == "REG;REG2"
    assert rs7["hop"] == 1
