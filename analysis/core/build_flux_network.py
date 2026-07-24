#!/usr/bin/env python3
"""Build one tissue's disease-regulatory flux network (node roles + flux edges).

Consumes the paper-exact flux edges written by bin/build_flux_map.py
(flux_<tissue>_<trait>.csv: regulator -> disease-gene attribution flux, with
sum_A flux(A -> G) = z_trans(G) exactly), so this network and the flux-edge
supplementary table are computed from one and the same flux. Nodes are all
disease-significant genes (association p_adj < FDR) together with every regulator
that contributes flux to a significant gene.

Per-gene flux quantities:
  in_degree  = # distinct regulators feeding the gene (== # incoming flux edges)
  in_flux    = sum |flux(A->gene)|                 (received, as a target)
  out_degree = # distinct sig targets the gene feeds
  out_flux   = sum |flux(gene->G)|                 (sent; ==0 unless the gene is a
               cis-eGene regulator -- only cis-eGenes can broadcast)
  net_in     = sum flux(A->gene)                   (signed; == z_trans of a sig gene)
  net_out    = sum flux(gene->G)                   (signed; for PMR coherence)
  influx_ratio = in_flux / (in_flux + out_flux)    (1 when out_flux==0; NaN if both 0)
  dominant_edge_frac = max_A |flux(A->gene)| / in_flux   (top regulator's share of the
               received flux; 1 when in_degree==1, NaN when in_flux==0)
  single_edge_driven = disease_sig AND in_degree>=1 AND dominant_edge_frac >= --single-edge-frac
               -- one regulator carries most of the trans signal (not distributed convergence).

Classification (robust to the structural out_flux==0):
  CORE       = disease-sig AND [ (out_flux==0 AND in_degree >= k) OR
                                 (out_flux>0  AND influx_ratio > 0.7) ]
  PERIPHERAL = out_flux>0 AND influx_ratio < 0.3
  MIXED      = out_flux>0 AND 0.3 <= influx_ratio <= 0.7
  PMR        = PERIPHERAL AND permutation-significant net out-flux (sign-flip null)
k defaults to the median in_degree among sig genes that receive (in_degree>=1).

Outputs (per tissue, per trait):
  flux_edges_<T>_<trait>.tsv   tissue regulator target flux abs_flux target_sig
  flux_nodes_<T>_<trait>.tsv   one row per gene: degrees/fluxes/ratio/class/PMR,
                               dominant_edge_frac and single_edge_driven
"""
import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

FDR = 0.05


def signflip_z_p(flux_vec):
    """Analytic sign-flip (Rademacher) null for net-flux coherence: under random signs the net
    has mean 0 and variance sum f_i^2, so z = net / sqrt(sum f_i^2) ~ N(0,1). A perfectly coherent
    regulator with k equal-magnitude targets gives z = sqrt(k)."""
    k = len(flux_vec)
    if k < 2:
        return np.nan, np.nan
    s2 = float(np.sum(flux_vec ** 2))
    if s2 == 0:
        return np.nan, np.nan
    z = float(flux_vec.sum()) / np.sqrt(s2)
    return z, float(2 * norm.sf(abs(z)))


def bh_q(pseries):
    """BH FDR over a Series of p-values (index preserved), NaNs passed through."""
    p = pseries.dropna().sort_values()
    m = len(p)
    if m == 0:
        return pd.Series(np.nan, index=pseries.index)
    q = p.values * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1].clip(max=1.0)
    qmap = dict(zip(p.index, q))
    return pseries.index.to_series().map(qmap)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--trait", required=True)
    ap.add_argument("--flux", required=True,
                    help="flux_<t>_<trait>.csv from bin/build_flux_map.py (the paper-exact edges)")
    ap.add_argument("--assoc", required=True, help="association_<t>_<trait>.tsv")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--core-k", default="median", help="'median' or an integer in-degree floor")
    ap.add_argument("--pmr-min-outdeg", type=int, default=5,
                    help="PMR candidate: min # disease targets a peripheral gene feeds")
    ap.add_argument("--pmr-min-coherence", type=float, default=0.6,
                    help="PMR candidate: min |net_out|/out_flux sign coherence")
    ap.add_argument("--single-edge-frac", type=float, default=0.5,
                    help="single_edge_driven flag: min dominant_edge_frac to call a sig gene's "
                         "trans signal single-edge driven rather than multi-regulator convergence")
    args = ap.parse_args()
    t, trait = args.tissue, args.trait
    os.makedirs(args.out_dir, exist_ok=True)

    meta = pd.read_csv(args.assoc, sep="\t").set_index("gene_id")
    if "win_config" not in meta.columns and "config" in meta.columns:
        meta["win_config"] = meta["config"]
    sig = set(meta.index[meta["p_adj"] < FDR])

    # paper-exact flux edges (regulator -> disease gene), keyed on Ensembl gene ids
    fx = pd.read_csv(args.flux)
    E = pd.DataFrame({
        "tissue": t,
        "regulator": fx["regulator_gene_id"],
        "target": fx["target_gene_id"],
        "flux": fx["flux"].astype(float),
    })
    E["abs_flux"] = E["flux"].abs()
    E["target_sig"] = E["target"].isin(sig)
    ep = f"{args.out_dir}/flux_edges_{t}_{trait}.tsv"
    E.to_csv(ep, sep="\t", index=False)

    # ---- per-gene in (as target) / out (as regulator) ----
    indeg = E.groupby("target").agg(in_degree=("regulator", "nunique"),
                                    in_flux=("abs_flux", "sum"),
                                    net_in=("flux", "sum"))
    outdeg = E.groupby("regulator").agg(out_degree=("target", "nunique"),
                                        out_flux=("abs_flux", "sum"),
                                        net_out=("flux", "sum"))
    idx = sorted(set(E["regulator"]) | set(E["target"]) | sig)
    nodes = pd.DataFrame(index=idx).join(indeg).join(outdeg)
    for c in ["in_degree", "in_flux", "net_in", "out_degree", "out_flux", "net_out"]:
        nodes[c] = nodes[c].fillna(0.0)
    nodes["disease_sig"] = nodes.index.isin(sig)
    den = nodes["in_flux"] + nodes["out_flux"]
    nodes["influx_ratio"] = np.where(den > 0, nodes["in_flux"] / den, np.nan)

    # concentration of received flux: carried by ONE regulator or spread across many?
    max_in = E.groupby("target")["abs_flux"].max().reindex(nodes.index)
    nodes["dominant_edge_frac"] = np.where(nodes["in_flux"] > 0,
                                           max_in.to_numpy() / nodes["in_flux"], np.nan)
    nodes["single_edge_driven"] = (nodes["disease_sig"] & (nodes["in_degree"] >= 1)
                                   & (nodes["dominant_edge_frac"] >= args.single_edge_frac)
                                   ).fillna(False)

    # in-degree floor k for core
    if str(args.core_k).lower() == "median":
        recv = nodes.loc[nodes["disease_sig"] & (nodes["in_degree"] >= 1), "in_degree"]
        k = float(recv.median()) if len(recv) else 1.0
    else:
        k = float(args.core_k)

    out0 = nodes["out_flux"] == 0
    ratio = nodes["influx_ratio"]
    is_core = nodes["disease_sig"] & ((out0 & (nodes["in_degree"] >= k)) |
                                      (~out0 & (ratio > 0.7)))
    is_periph = (~out0) & (ratio < 0.3)
    is_mixed = (~out0) & (ratio >= 0.3) & (ratio <= 0.7)
    nodes["flux_class"] = np.select([is_core, is_periph, is_mixed],
                                    ["core", "peripheral", "mixed"], default="none")

    # ---- PMR: a peripheral broadcaster delivering COHERENT sign to MANY targets ----
    nodes["coherence"] = np.where(nodes["out_flux"] > 0,
                                  nodes["net_out"].abs() / nodes["out_flux"], np.nan)
    nodes["pmr_z"] = np.nan
    nodes["pmr_p"] = np.nan
    fl_by_reg = {r: grp["flux"].to_numpy() for r, grp in E.groupby("regulator")}
    for r in nodes.index[is_periph]:
        z, p = signflip_z_p(fl_by_reg[r])
        nodes.at[r, "pmr_z"] = z
        nodes.at[r, "pmr_p"] = p
    pmr_cand = (is_periph & (nodes["out_degree"] >= args.pmr_min_outdeg)
                & (nodes["coherence"] >= args.pmr_min_coherence))
    nodes["pmr_cand"] = pmr_cand
    nodes["pmr_q"] = bh_q(nodes.loc[pmr_cand, "pmr_p"]).reindex(nodes.index)
    nodes["is_pmr"] = pmr_cand & (nodes["pmr_p"] < 0.05)

    # attach association fields
    for col in ["z_cis", "z_trans", "z_twas", "p_adj", "win_config"]:
        if col in meta.columns:
            nodes[col] = nodes.index.map(meta[col].to_dict())

    nodes = nodes.reset_index().rename(columns={"index": "gene_id"})
    nodes.insert(0, "tissue", t)
    npth = f"{args.out_dir}/flux_nodes_{t}_{trait}.tsv"
    nodes.to_csv(npth, sep="\t", index=False)

    vc = nodes["flux_class"].value_counts().to_dict()
    n_sed = int(nodes["single_edge_driven"].sum())
    n_trans_hit = int((nodes["disease_sig"] & (nodes["in_degree"] >= 1)).sum())
    # exact-decomposition self-check: net_in of every sig gene equals its z_trans
    chk = nodes[nodes["disease_sig"] & (nodes["in_degree"] >= 1)]
    resid = float((chk["net_in"] - chk["z_trans"].astype(float)).abs().max()) if len(chk) else 0.0
    print(f"[{t}] edges={len(E)} nodes={len(nodes)} sig={len(sig)} k={k:.0f} | "
          f"core={vc.get('core',0)} peripheral={vc.get('peripheral',0)} "
          f"mixed={vc.get('mixed',0)} PMR={int(nodes['is_pmr'].sum())} "
          f"single_edge_driven={n_sed}/{n_trans_hit} | max|net_in-z_trans|={resid:.1e} "
          f"-> {ep} , {npth}")


if __name__ == "__main__":
    main()
