#!/usr/bin/env python3
"""Assemble the release supplementary tables from the per-tissue result files.

Each supplementary table combines the seven STARNET tissues into one file, with an
HGNC ``gene`` symbol column added and the gene_id kept as the disambiguator. The set
spans the expression models (S1), the CAD transcriptome-wide association (S2), the
regulatory network and flux map (S3-S5), replication and functional evidence
(S6-S8), and the per-SNP predictor weights (S9):

  Supplementary Table S1 -- per-gene channel-aware expression-model summary
      (selected configuration, per-channel LOO-R2 and ELPD, predictability).
      Source: results_cv/horseshoe_alltargets/models/<T>.horseshoe_metrics.tsv.gz

  Supplementary Table S2 -- transcriptome-wide association statistics for CAD
      (cis/trans decomposition of z_TWAS, p and FDR).
      Source: results_cv/horseshoe_alltargets/association/association_<T>_<trait>.tsv

  Supplementary Table S9 -- per-SNP weights of each gene's selected predictor,
      annotated with dbSNP rsID and GRCh37 coordinates/alleles so external cohorts
      can apply them, and, for trans SNPs, the regulator gene(s) the SNP is a
      cis-eQTL of. Built only when --variant-map is given.
      Source: results_cv/horseshoe_alltargets/models/<T>.horseshoe_weights.tsv.gz

The tables carry gene-, tissue-, edge- and per-SNP-weight summaries only; no
individual genotypes, expression values or sample identifiers are included, so none
redistributes restricted STARNET data. S9 identifies each SNP by its public dbSNP
rsID and alleles.

A README data dictionary describing every column is written alongside them.

Input result locations resolve under RESULTS (default results_cv/, override with FLUX_RESULTS);
gene symbols come from the annotation in $GENE_ANNOT (or --gene-annot).

Usage:
  GENE_ANNOT=<gencode.v19.genes.tsv> python analysis/tables/build_supp_tables.py \
      [--trait CAD_aragam2022] [--outdir supplementary_tables] \
      [--gene-annot <gencode.v19.genes.tsv>] \
      [--variant-map <variant annotation TSV>] [--variant-map-id-col original_id]
"""
import argparse
import gzip
import os
import re
import sys
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_ROOT / "utils"))
import gene_labels as gl

TISSUES = ["AOR", "Blood", "LIV", "MAM", "SF", "SKLM", "VAF"]
# Input result locations resolve under RESULTS (default results_cv/, override with FLUX_RESULTS).
RESULTS = os.environ.get("FLUX_RESULTS", str(_ROOT / "results_cv"))
MODELS = f"{RESULTS}/horseshoe_alltargets/models"
WEIGHTS = f"{RESULTS}/horseshoe_alltargets/models/{{t}}.horseshoe_weights.tsv.gz"
SNP_SD = f"{RESULTS}/horseshoe_alltargets/models/{{t}}.snp_sd.tsv.gz"
ASSOC = f"{RESULTS}/horseshoe_alltargets/association"
GRN = f"{RESULTS}/network_alltargets/dag_{{t}}_orig_kde_fdr15_alltargets.csv"
FLUX = f"{RESULTS}/horseshoe_alltargets/flux/flux_edges_{{t}}_{{trait}}.tsv"
FEATURES = f"{RESULTS}/features_alltargets/{{t}}.trans_features.tsv.gz"
NODES = f"{RESULTS}/horseshoe_alltargets/flux/flux_nodes_{{t}}_{{trait}}.tsv"
FINNGEN = f"{RESULTS}/horseshoe_alltargets/omnipath_validation/finngen_replication.tsv"
FUNCENR = f"{RESULTS}/horseshoe_alltargets/omnipath_validation/flux_function_enrichment.tsv"
# Curated gene-set libraries (.gmt) ship in the repo under data/genesets/.
GENESETS = str(_ROOT / "data" / "genesets")
GMT_LIBS = [
    ("MGI", "MGI_Mammalian_Phenotype_Level_4_2021.gmt"),
    ("GOBP", "GO_Biological_Process_2023.gmt"),
    ("REAC", "Reactome_2022.gmt"),
    ("KEGG", "KEGG_2021_Human.gmt"),
    ("DIS", "Jensen_DISEASES.gmt"),
]
CV_KEYWORDS = re.compile(
    r"cardiovascular|cardiac|heart|myocard|cardiomyo|vascul|vasocon|vasodil|angiogen|"
    r"endotheli|blood vessel|artery|arterial|aorta|aortic|atheroscler|smooth muscle|"
    r"blood pressure|hypertens|coagul|platelet|thrombo|hemostasis|fibrin|lipid|lipoprotein|"
    r"cholesterol|choleste|triglycer|sterol|foam cell|coronary|ischemi|infarct|elastic fiber|"
    r"extracellular matrix", re.I)

# Release-safe column selection, in output order (gene_id kept as disambiguator).
S1_COLS = [
    "gene", "gene_id", "tissue", "gene_class",
    "n_cis", "n_trans",
    "loo_r2_cis", "loo_r2_trans", "loo_r2_cistrans",
    "elpd_cis", "elpd_trans", "elpd_cistrans", "delta_elpd",
    "best_config", "best_loo_r2", "predictable",
]
S2_COLS = [
    "gene", "gene_id", "tissue", "win_method", "win_config",
    "n_cis", "n_trans",
    "z_cis", "z_trans", "z_twas",
    "sigma_g", "p", "p_adj",
]

S1_DICT = [
    ("gene", "HGNC gene symbol (GENCODE v19); the base ENSG is shown when no symbol maps"),
    ("gene_id", "Ensembl gene identifier (versioned ENSG), the unambiguous key"),
    ("tissue", "Tissue: AOR, Blood, LIV, MAM, SF, SKLM, VAF"),
    ("gene_class", "Channels available for the gene: cis_only, trans_only or cis_trans"),
    ("n_cis", "Number of cis-eQTL SNPs in the cis channel"),
    ("n_trans", "Number of trans-feature SNPs (regulators' cis-eQTLs) in the trans channel"),
    ("loo_r2_cis", "Leave-one-out R2 of the cis-only model (blank if unavailable)"),
    ("loo_r2_trans", "Leave-one-out R2 of the trans-only model (blank if unavailable)"),
    ("loo_r2_cistrans", "Leave-one-out R2 of the combined cis+trans model (blank if unavailable)"),
    ("elpd_cis", "PSIS-LOO expected log pointwise predictive density, cis-only model"),
    ("elpd_trans", "PSIS-LOO ELPD, trans-only model"),
    ("elpd_cistrans", "PSIS-LOO ELPD, combined cis+trans model"),
    ("delta_elpd", "ELPD(cis+trans) - ELPD(cis-only); positive favours the combined model"),
    ("best_config", "Selected configuration by ELPD: cis_only, trans_only or cis_trans"),
    ("best_loo_r2", "Leave-one-out R2 of the selected configuration"),
    ("predictable", "True when best_loo_r2 >= 0.01 (used for downstream association)"),
]
S3_COLS = [
    "source", "source_gene_id", "target", "target_gene_id", "tissue",
    "probability", "qvalue",
]

S3_DICT = [
    ("source", "HGNC symbol of the regulator (source) gene"),
    ("source_gene_id", "Ensembl gene identifier of the regulator"),
    ("target", "HGNC symbol of the target gene"),
    ("target_gene_id", "Ensembl gene identifier of the target"),
    ("tissue", "Tissue: AOR, Blood, LIV, MAM, SF, SKLM, VAF"),
    ("probability", "Posterior probability of the directed edge (Findr edge score)"),
    ("qvalue", "Global false-discovery-rate q-value for the edge"),
]

S4_COLS = [
    "regulator", "regulator_gene_id", "target", "target_gene_id", "tissue", "hop", "flux",
]

S4_DICT = [
    ("regulator", "HGNC symbol of the upstream regulator gene"),
    ("regulator_gene_id", "Ensembl gene identifier of the regulator"),
    ("target", "HGNC symbol of the target gene (every target is disease-significant, FDR < 0.05)"),
    ("target_gene_id", "Ensembl gene identifier of the target"),
    ("tissue", "Tissue: AOR, Blood, LIV, MAM, SF, SKLM, VAF"),
    ("hop", "Network distance from regulator to target: 1 = direct parent, 2 = grandparent"),
    ("flux", "Signed disease flux the regulator delivers to the target; positive drives the "
             "target's expression toward risk-increasing, and the fluxes into a target sum to "
             "its trans association component"),
]

S5_COLS = [
    "gene", "gene_id", "tissue", "role", "disease_sig", "z_twas",
    "in_degree", "out_degree", "n_disease_targets", "in_flux", "out_flux", "is_pmr",
]

S5_DICT = [
    ("gene", "HGNC gene symbol (GENCODE v19); the base ENSG is shown when no symbol maps"),
    ("gene_id", "Ensembl gene identifier (versioned ENSG), the unambiguous key"),
    ("tissue", "Tissue: AOR, Blood, LIV, MAM, SF, SKLM, VAF"),
    ("role", "Node role in this tissue: 'core master regulator' (disease-significant and "
             "broadcasts to > 10 disease genes), 'sink' (disease-significant, not a broad "
             "broadcaster), 'key driver' (not significant but broadcasts to > 10 disease genes), "
             "or 'minor regulator' (neither)"),
    ("disease_sig", "True when the gene is disease-significant (FDR < 0.05) in this tissue"),
    ("z_twas", "TWAS association statistic for the gene in this tissue"),
    ("in_degree", "Number of upstream regulators that feed the gene"),
    ("out_degree", "Number of genes the gene feeds"),
    ("n_disease_targets", "Number of disease-significant genes the gene broadcasts to (breadth "
                          "used to assign the role)"),
    ("in_flux", "Total absolute flux the gene receives from its regulators"),
    ("out_flux", "Total absolute flux the gene broadcasts to disease genes"),
    ("is_pmr", "True when the gene is a peripheral master regulator (drives many disease genes "
               "coherently beyond chance)"),
]

S8_COLS = ["gene_set", "lib", "term", "p", "adj_p", "combined", "n_overlap"]

S8_DICT = [
    ("gene_set", "Gene set tested: cis_sig (significant genes won by a cis-only or cis+trans "
                 "predictor), trans_only (significant genes won by a trans-only predictor), or "
                 "contributor (non-significant regulators broadcasting flux to significant genes)"),
    ("lib", "Annotation library: GO Biological Process 2023 or Reactome 2022"),
    ("term", "Gene-set term name"),
    ("p", "Enrichment p-value (Fisher exact test)"),
    ("adj_p", "Benjamini-Hochberg adjusted p-value across terms (enriched at < 0.05)"),
    ("combined", "Enrichr combined score"),
    ("n_overlap", "Number of genes in the set that overlap the term"),
]

S7_COLS = ["gene", "gene_id", "in_cv_universe", "in_mgi_ko_cv"]

S7_DICT = [
    ("gene", "HGNC gene symbol; one row per tested gene"),
    ("gene_id", "Representative Ensembl gene identifier (version-stripped)"),
    ("in_cv_universe", "True if the gene is annotated to any cardiovascular term across the MGI, "
                       "GO Biological Process, Reactome, KEGG or Jensen DISEASES libraries "
                       "(the curated cardiovascular universe)"),
    ("in_mgi_ko_cv", "True if the gene is annotated to an MGI mouse-knockout cardiovascular "
                     "phenotype (a readout independent of human genetics)"),
]

S6_COLS = [
    "gene", "gene_id", "tissue", "win_config", "z_discovery", "z_finngen", "concordant",
]

S6_DICT = [
    ("gene", "HGNC gene symbol (GENCODE v19); the base ENSG is shown when no symbol maps"),
    ("gene_id", "Ensembl gene identifier (versioned ENSG), the unambiguous key"),
    ("tissue", "Tissue: AOR, Blood, LIV, MAM, SF, SKLM, VAF"),
    ("win_config", "Selected predictor configuration: cis_only, trans_only or cis_trans"),
    ("z_discovery", "Discovery TWAS statistic (CARDIoGRAMplusC4D CAD GWAS)"),
    ("z_finngen", "TWAS statistic recomputed from FinnGen R11 coronary heart disease (I9_CHD) "
                  "using the same predictor weights"),
    ("concordant", "True when the FinnGen direction matches the discovery direction "
                   "(sign of z_finngen equals sign of z_discovery)"),
]

S2_DICT = [
    ("gene", "HGNC gene symbol (GENCODE v19); the base ENSG is shown when no symbol maps"),
    ("gene_id", "Ensembl gene identifier (versioned ENSG), the unambiguous key"),
    ("tissue", "Tissue: AOR, Blood, LIV, MAM, SF, SKLM, VAF"),
    ("win_method", "Expression-model method of the selected predictor (horseshoe)"),
    ("win_config", "Selected predictor configuration: cis_only, trans_only or cis_trans"),
    ("n_cis", "Number of cis SNPs contributing to the statistic"),
    ("n_trans", "Number of trans SNPs contributing to the statistic"),
    ("z_cis", "Cis component of the TWAS statistic"),
    ("z_trans", "Trans component of the TWAS statistic (z_twas = z_cis + z_trans)"),
    ("z_twas", "Summary-statistic TWAS z against the CAD GWAS"),
    ("sigma_g", "Standard deviation of genetically predicted expression (statistic denominator)"),
    ("p", "Two-sided normal p-value for z_twas"),
    ("p_adj", "Benjamini-Hochberg FDR-adjusted p-value (per tissue); significant at < 0.05"),
]

S9_COLS = [
    "gene", "gene_id", "tissue", "config", "channel", "rs_id",
    "chromosome", "position", "ref", "alt", "weight", "sd",
    "source_gene", "source_gene_id", "hop",
]

S9_DICT = [
    ("gene", "HGNC gene symbol (GENCODE v19); the base ENSG is shown when no symbol maps"),
    ("gene_id", "Ensembl gene identifier (versioned ENSG), the unambiguous key"),
    ("tissue", "Tissue: AOR, Blood, LIV, MAM, SF, SKLM, VAF"),
    ("config", "Selected predictor configuration the weight belongs to: cis_only, "
               "trans_only or cis_trans"),
    ("channel", "Channel of the SNP within the predictor: cis (the gene's own cis-eQTL) "
                "or trans (a regulator's cis-eQTL)"),
    ("rs_id", "dbSNP rsID of the predictor SNP"),
    ("chromosome", "Chromosome (GRCh37)"),
    ("position", "Base-pair position (GRCh37)"),
    ("ref", "Reference allele"),
    ("alt", "Alternative allele; the weight is applied to the ALT-allele dosage"),
    ("weight", "Predictor weight multiplying the ALT-allele dosage of this SNP; the "
               "genetically predicted expression is the sum of weight x dosage over the "
               "gene's SNPs (one row per SNP)"),
    ("sd", "Standard deviation of the SNP's ALT-dosage over the training samples; used for "
           "summary-statistic TWAS (z_TWAS = sum_j weight_j x sd_j x z_gwas_j / sigma_g, with "
           "sigma_g from Supplementary Table S2). Blank for the few SNPs outside the "
           "genotype panel used for the association"),
    ("source_gene", "For a trans SNP, the HGNC symbol(s) of the regulator gene(s) the SNP "
                    "is a cis-eQTL of; ';'-separated when the SNP tags several regulators. "
                    "Blank for cis SNPs"),
    ("source_gene_id", "Ensembl gene identifier(s) of the regulator(s), matching source_gene"),
    ("hop", "Network distance from the regulator to the gene for a trans SNP (1 = direct "
            "parent, 2 = grandparent; the minimum when several regulators are tagged); "
            "blank for cis SNPs"),
]


def load_models(mapping):
    """Combine the seven per-tissue expression-model metrics into one table."""
    frames = []
    for t in TISSUES:
        path = f"{MODELS}/{t}.horseshoe_metrics.tsv.gz"
        df = pd.read_csv(path, sep="\t")
        df.insert(0, "tissue", t)
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out = gl.add_symbol_column(out, mapping, id_col="gene_id", out_col="gene")
    return out[S1_COLS]


def load_grn(mapping):
    """Combine the seven per-tissue reconstructed GRN edge lists into one table."""
    frames = []
    for t in TISSUES:
        df = pd.read_csv(GRN.format(t=t), usecols=["Source", "Target", "Probability", "qvalue"])
        df = df.rename(columns={"Source": "source_gene_id", "Target": "target_gene_id",
                                "Probability": "probability"})
        df.insert(0, "tissue", t)
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out["source"] = out["source_gene_id"].map(lambda g: gl.resolve(str(g), mapping))
    out["target"] = out["target_gene_id"].map(lambda g: gl.resolve(str(g), mapping))
    return out[S3_COLS]


def load_flux(mapping, trait):
    """Combine the seven per-tissue flux-edge tables into one table, labelling each
    regulator-to-target edge with its network distance (hop 1 = parent, 2 = grandparent)."""
    frames = []
    for t in TISSUES:
        df = pd.read_csv(FLUX.format(t=t, trait=trait), sep="\t",
                         usecols=["tissue", "regulator", "target", "flux"])
        tf = pd.read_csv(FEATURES.format(t=t), sep="\t",
                         usecols=["gene_id", "hop", "source_gene_id"]).drop_duplicates()
        hop = tf.groupby(["gene_id", "source_gene_id"])["hop"].min().reset_index()
        df = df.merge(hop, left_on=["target", "regulator"],
                      right_on=["gene_id", "source_gene_id"], how="left")
        df = df.rename(columns={"regulator": "regulator_gene_id", "target": "target_gene_id"})
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out["hop"] = out["hop"].astype("Int64")
    out["regulator"] = out["regulator_gene_id"].map(lambda g: gl.resolve(str(g), mapping))
    out["target"] = out["target_gene_id"].map(lambda g: gl.resolve(str(g), mapping))
    return out[S4_COLS]


def load_functional_enrichment():
    """Load the per-term GO/Reactome enrichment of the flux-map gene sets."""
    df = pd.read_csv(FUNCENR, sep="\t").rename(columns={"set": "gene_set"})
    return df[S8_COLS]


def _load_gmt(path, prefix):
    """Read an Enrichr .gmt gene-set library: term (prefixed) -> set of upper-case symbols."""
    d = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        d[f"{prefix}:{f[0]}"] = set(g.upper() for g in f[2:] if g)
    return d


def load_cv_universe(mapping, trait):
    """Build the curated cardiovascular gene universe and the MGI mouse-knockout subset,
    restricted to the genes tested across the seven tissues, as per-gene membership flags."""
    ids = set()
    for t in TISSUES:
        a = pd.read_csv(f"{ASSOC}/association_{t}_{trait}.tsv", sep="\t", usecols=["gene_id"])
        ids.update(str(g).split(".")[0] for g in a["gene_id"])
    sym_of = {g: gl.resolve(g, mapping) for g in ids}
    tested_sym = set(sym_of.values())

    cvterms = {}
    for prefix, fn in GMT_LIBS:
        for term, genes in _load_gmt(f"{GENESETS}/{fn}", prefix).items():
            if CV_KEYWORDS.search(term):
                cvterms[term] = genes
    cv = set().union(*cvterms.values()) & tested_sym
    mgi = set().union(*[g for term, g in cvterms.items() if term.startswith("MGI")]) & tested_sym

    rep = {}
    for g in sorted(ids):
        rep.setdefault(sym_of[g], g)
    rows = [{"gene": s, "gene_id": rep[s],
             "in_cv_universe": s in cv, "in_mgi_ko_cv": s in mgi}
            for s in sorted(tested_sym)]
    return pd.DataFrame(rows)[S7_COLS]


def load_finngen(mapping):
    """Load the per-gene FinnGen R11 CHD replication (discovery vs FinnGen TWAS direction)."""
    df = pd.read_csv(FINNGEN, sep="\t")
    df = df.rename(columns={"wc": "win_config", "zd": "z_discovery", "zfin": "z_finngen"})
    df["concordant"] = df["concord"] == 1
    df = gl.add_symbol_column(df, mapping, id_col="gene_id", out_col="gene")
    return df[S6_COLS]


def load_node_roles(mapping, trait):
    """Combine the seven per-tissue flux-node tables and assign each gene its four-class role
    in that tissue (core master regulator / sink / key driver / minor regulator)."""
    frames = []
    for t in TISSUES:
        frames.append(pd.read_csv(NODES.format(t=t, trait=trait), sep="\t"))
    n = pd.concat(frames, ignore_index=True)
    n["n_disease_targets"] = pd.to_numeric(n["n_disease_targets"], errors="coerce").fillna(0)
    n["disease_sig"] = n["disease_sig"] == True
    sig, bcast = n["disease_sig"], n["n_disease_targets"] > 10
    role = pd.Series("minor regulator", index=n.index)
    role[sig & ~bcast] = "sink"
    role[~sig & bcast] = "key driver"
    role[sig & bcast] = "core master regulator"
    n["role"] = role
    n["gene"] = n["gene_id"].map(lambda g: gl.resolve(str(g), mapping))
    return n[S5_COLS]


def load_assoc(mapping, trait):
    """Combine the seven per-tissue CAD association tables into one table."""
    frames = []
    for t in TISSUES:
        path = f"{ASSOC}/association_{t}_{trait}.tsv"
        frames.append(pd.read_csv(path, sep="\t"))
    out = pd.concat(frames, ignore_index=True)
    out = gl.add_symbol_column(out, mapping, id_col="gene_id", out_col="gene")
    return out[S2_COLS]


def load_snp_coords(path, id_col, wanted):
    """Look up chromosome/position/ref/alt for the requested variant IDs from a variant
    annotation file, streaming so only the requested IDs are held in memory. Returns a
    dict id -> (chromosome, position, ref, alt)."""
    opener = gzip.open if path.endswith(".gz") else open
    coords = {}
    with opener(path, "rt") as fh:
        idx = {c: i for i, c in enumerate(fh.readline().rstrip("\n").split("\t"))}
        i_id = idx[id_col]
        i_chr, i_pos, i_ref, i_alt = idx["chromosome"], idx["position"], idx["ref"], idx["alt"]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[i_id] in wanted:
                coords[f[i_id]] = (f[i_chr], f[i_pos], f[i_ref], f[i_alt])
    return coords


def load_weights(mapping, variant_map, id_col):
    """Combine the per-tissue predictor weights into one table, keeping only the weights
    of each gene's selected configuration for predictable genes, annotating every SNP with
    its dbSNP coordinates and alleles from the variant map, and, for trans SNPs, the
    regulator gene(s) the SNP is a cis-eQTL of (from the trans-feature tables)."""
    frames = []
    for t in TISSUES:
        w = pd.read_csv(WEIGHTS.format(t=t), sep="\t")
        w = w[w["config"].isin(["cis_only", "trans_only", "cis_trans"])]
        m = pd.read_csv(f"{MODELS}/{t}.horseshoe_metrics.tsv.gz", sep="\t",
                        usecols=["gene_id", "best_config", "predictable"])
        m = m[m["predictable"].astype(str).str.lower() == "true"]
        sel = w.merge(m[["gene_id", "best_config"]], on="gene_id")
        sel = sel[sel["config"] == sel["best_config"]]
        sel.insert(0, "tissue", t)

        # For trans SNPs, attach the regulator(s) whose cis-eQTL the SNP is. A SNP can tag
        # several regulators (and at different hops), so collapse them to one row per SNP:
        # ';'-joined regulator ids and the minimum hop.
        tf = pd.read_csv(FEATURES.format(t=t), sep="\t",
                         usecols=["gene_id", "hop", "source_gene_id", "rs_id"])
        reg = tf.groupby(["gene_id", "rs_id"]).agg(
            source_gene_id=("source_gene_id", lambda s: ";".join(sorted(set(s)))),
            hop=("hop", "min")).reset_index()
        sel = sel.merge(reg, on=["gene_id", "rs_id"], how="left")
        cis = sel["channel"] != "trans"
        sel.loc[cis, ["source_gene_id", "hop"]] = pd.NA

        # Per-SNP ALT-dosage SD over the training samples (the summary-TWAS variance term).
        sdf = pd.read_csv(SNP_SD.format(t=t), sep="\t", usecols=["rs_id", "sd"])
        sel = sel.merge(sdf, on="rs_id", how="left")

        frames.append(sel[["gene_id", "tissue", "config", "channel", "rs_id",
                            "weight", "sd", "source_gene_id", "hop"]])
    out = pd.concat(frames, ignore_index=True)
    out["hop"] = out["hop"].astype("Int64")

    coords = load_snp_coords(variant_map, id_col, set(out["rs_id"]))
    missing = sorted(set(out["rs_id"]) - set(coords))
    if missing:
        print(f"S9 warning: {len(missing)} SNP(s) absent from the variant map "
              f"(alleles left blank), e.g. {missing[:3]}", file=sys.stderr)
    for i, col in enumerate(["chromosome", "position", "ref", "alt"]):
        out[col] = out["rs_id"].map(lambda v, i=i: coords.get(v, (None, None, None, None))[i])
    out["source_gene"] = out["source_gene_id"].map(
        lambda s: "" if pd.isna(s) else ";".join(gl.resolve(g, mapping) for g in s.split(";")))
    out = gl.add_symbol_column(out, mapping, id_col="gene_id", out_col="gene")
    return out[S9_COLS]


def write_readme(path, trait, n1, n2, n3, n4, n5, n6, n7, n8, n9=None):
    """Write the data-dictionary README for the tables."""
    def block(title, dictionary):
        lines = [title, "", "Columns:"]
        lines += [f"  {name}\t{desc}" for name, desc in dictionary]
        return "\n".join(lines)

    sections = [
        "Supplementary tables for the disease-regulatory flux map (CAD)",
        "Each table combines the seven tissues (AOR, Blood, LIV, MAM, SF, "
        "SKLM, VAF). The tables carry gene-, edge- and tissue-level summaries "
        "only and contain no individual-level genotype or expression data.",
        "Network vs flux: S3 is the complete regulatory network (every directed "
        "gene-to-gene edge), from which the up-to-two-hop ancestry of any gene is "
        "recovered (parents = direct edges into the gene; grandparents = parents of "
        "those parents). S4 is the disease-specific flux on top of it: each ancestor "
        "is attributed directly to the disease gene it feeds (hop 1 or 2), so the "
        "grandparent-to-parent routing edge is read from S3 while its signed flux is "
        "read from S4.",
        block(
            f"Supplementary Table S1 -- per-gene channel-aware expression-model "
            f"summary ({n1} gene-tissue rows).",
            S1_DICT,
        ),
        block(
            f"Supplementary Table S2 -- CAD transcriptome-wide association "
            f"statistics ({n2} gene-tissue rows).",
            S2_DICT,
        ),
        block(
            f"Supplementary Table S3 -- reconstructed gene regulatory network "
            f"edges ({n3} directed edges).",
            S3_DICT,
        ),
        block(
            f"Supplementary Table S4 -- disease-regulatory flux edges "
            f"({n4} regulator-to-target edges).",
            S4_DICT,
        ),
        block(
            f"Supplementary Table S5 -- per-gene-tissue flux node roles "
            f"({n5} gene-tissue rows).",
            S5_DICT,
        ),
        block(
            f"Supplementary Table S6 -- FinnGen R11 CHD replication of the "
            f"associations ({n6} gene-tissue rows).",
            S6_DICT,
        ),
        block(
            f"Supplementary Table S7 -- cardiovascular gene-set membership "
            f"({n7} tested genes).",
            S7_DICT,
        ),
        block(
            f"Supplementary Table S8 -- per-term functional enrichment of the "
            f"flux-map gene sets ({n8} term tests).",
            S8_DICT,
        ),
    ]
    if n9 is not None:
        sections.append(block(
            f"Supplementary Table S9 -- per-SNP weights of each gene's selected "
            f"predictor ({n9} SNP-gene-tissue rows).",
            S9_DICT,
        ))
    text = "\n\n".join(sections) + "\n"
    with open(path, "w") as fh:
        fh.write(text)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trait", default="CAD_aragam2022")
    ap.add_argument("--outdir", default="supplementary_tables")
    ap.add_argument("--gene-annot", default=None,
                    help="GENCODE v19 gene annotation TSV for HGNC symbols "
                         "(default: the $GENE_ANNOT environment variable).")
    ap.add_argument("--variant-map",
                    help="TSV annotating the predictor variant IDs with chromosome/position/"
                         "ref/alt. When given, Supplementary Table S9 (per-SNP predictor "
                         "weights with alleles) is built.")
    ap.add_argument("--variant-map-id-col", default="original_id",
                    help="Column in --variant-map holding the same variant IDs as the "
                         "weights' rs_id column (default: original_id).")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    mapping = gl.load_label_map(gl.default_annot(args.gene_annot))

    s1 = load_models(mapping)
    s2 = load_assoc(mapping, args.trait)
    s3 = load_grn(mapping)
    s4 = load_flux(mapping, args.trait)
    s5 = load_node_roles(mapping, args.trait)
    s6 = load_finngen(mapping)
    s7 = load_cv_universe(mapping, args.trait)
    s8 = load_functional_enrichment()
    s9 = load_weights(mapping, args.variant_map, args.variant_map_id_col) \
        if args.variant_map else None

    p1 = f"{args.outdir}/Supplementary_Table_S1_expression_models.csv"
    p2 = f"{args.outdir}/Supplementary_Table_S2_twas_association.csv"
    p3 = f"{args.outdir}/Supplementary_Table_S3_grn_edges.csv"
    p4 = f"{args.outdir}/Supplementary_Table_S4_flux_edges.csv"
    p5 = f"{args.outdir}/Supplementary_Table_S5_flux_node_roles.csv"
    p6 = f"{args.outdir}/Supplementary_Table_S6_finngen_replication.csv"
    p7 = f"{args.outdir}/Supplementary_Table_S7_cardiovascular_geneset.csv"
    p8 = f"{args.outdir}/Supplementary_Table_S8_functional_enrichment.csv"
    s1.to_csv(p1, index=False)
    s2.to_csv(p2, index=False)
    s3.to_csv(p3, index=False)
    s4.to_csv(p4, index=False)
    s5.to_csv(p5, index=False)
    s6.to_csv(p6, index=False)
    s7.to_csv(p7, index=False)
    s8.to_csv(p8, index=False)
    if s9 is not None:
        p9 = f"{args.outdir}/Supplementary_Table_S9_predictor_weights.csv.gz"
        s9.to_csv(p9, index=False, compression="gzip")
    write_readme(f"{args.outdir}/README.md", args.trait,
                 len(s1), len(s2), len(s3), len(s4), len(s5), len(s6), len(s7), len(s8),
                 None if s9 is None else len(s9))

    print(f"S1 expression models: {len(s1)} rows -> {p1}")
    print(f"S2 TWAS association:  {len(s2)} rows -> {p2}")
    print(f"S3 GRN edges:         {len(s3)} rows -> {p3}")
    print(f"S4 flux edges:        {len(s4)} rows -> {p4}")
    print(f"S5 flux node roles:   {len(s5)} rows -> {p5}")
    print(f"S6 FinnGen replication: {len(s6)} rows -> {p6}")
    print(f"S7 cardiovascular geneset: {len(s7)} rows -> {p7}")
    print(f"S8 functional enrichment: {len(s8)} rows -> {p8}")
    if s9 is not None:
        print(f"S9 predictor weights: {len(s9)} rows -> {p9}")
    else:
        print("S9 predictor weights: skipped (no --variant-map given)")
    print(f"README data dictionary -> {args.outdir}/README.md")


if __name__ == "__main__":
    main()
