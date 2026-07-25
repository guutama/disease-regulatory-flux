#!/usr/bin/env python3
"""Interactive Pyvis (vis.js) rendering of the largest connected component of a tissue's flux map.

One self-contained HTML per tissue. Disease genes (purple, sized by |z_trans|, labelled with
chromosome + cis/trans TWAS Z), key drivers -- non-significant regulators
that feed >= --min-targets disease genes (teal, sized by how many they feed) -- and private
regulators (light blue); a GOLD RING marks a core master regulator (a significant gene that itself
feeds >= --min-targets disease genes). Each
node label ends with its chromosome so same-chromosome regulators are easy to spot. Edges carry the
signed flux (red toward risk, blue away), width ~ |flux|, directed; an edge is dashed when the
regulator is on the SAME chromosome as its target (genomic distance in Mb shown -- a co-expression
clue). Full-viewport, draggable/zoomable, with PDF/PNG download buttons.

The displayed component is the largest connected component of the whole flux map -- every disease
gene together with its up-to-2-hop GRN ancestors (regulators) and all the flux edges among them.

Inputs (per tissue key; keys discovered from --flux-dir unless --keys is given):
  --flux-dir   flux_nodes_<key>_<trait>.tsv  and  flux_edges_<key>_<trait>.tsv
  --gencode    GENCODE genes tsv (gene_id, gene_name, chr, start, end) for symbols + chromosome

Output (--outdir): fig_flux_lcc_<key>.html (one per tissue).

Usage:
  plot_flux_lcc_network.py --flux-dir DIR --gencode FILE --outdir DIR \
      [--trait CAD_aragam2022] [--keys AOR LIV ...] [--conv-min 0.7] [--min-in-degree 2]
"""
import argparse
import glob
import os
import sys
from pathlib import Path

import pandas as pd
import networkx as nx
from pyvis.network import Network

_REPO = Path(__file__).resolve().parents[2]                  # analysis/figure_scripts/ -> repo root
sys.path.insert(0, str(_REPO / "utils"))
import gene_labels as gl

C_CORE, C_POS, C_NEG, C_KD, C_PRIV, C_DG = "#6a51a3", "#c1121f", "#2c6fbb", "#0d9488", "#9ecae1", "#e8a000"
# nicer per-tissue display names (fall back to the key itself)
TISSUE_NAME = {"AOR": "Aorta", "MAM": "Mammary artery", "LIV": "Liver", "VAF": "Visceral fat",
               "SF": "Subcutaneous fat", "SKLM": "Skeletal muscle", "Blood": "Blood"}

LEGEND = """
<style>html,body{margin:0;padding:0;overflow:hidden;} .card,.container,#mynetwork{margin:0!important;padding:0!important;border:0!important;}
#mynetwork{width:100vw!important;height:100vh!important;}</style>
<div style="position:absolute;top:12px;left:12px;background:rgba(255,255,255,0.94);border:1px solid #ccc;
 border-radius:8px;padding:12px 16px;font-family:Arial,sans-serif;font-size:13px;z-index:999;line-height:1.7;">
 <div style="font-size:15px;font-weight:bold;margin-bottom:4px;">LABEL flux network</div>
 <div style="color:#666;margin-bottom:6px;">largest connected component: NDG disease genes, NE edges &middot; labels end with chromosome</div>
 <span style="color:#6a51a3;font-size:22px;vertical-align:middle;">&#9679;</span> disease gene (size ~ |Z<sub>trans</sub>|; shows chr, Z<sub>cis</sub>/Z<sub>trans</sub>)<br>
 <span style="color:#0d9488;font-size:17px;vertical-align:middle;">&#9679;</span> key driver: non-significant regulator feeding &ge;MINTGT disease genes (size ~ # targets)<br>
 <span style="color:#9ecae1;font-size:14px;vertical-align:middle;">&#9679;</span> private regulator (feeds &lt;MINTGT disease genes)<br>
 <span style="display:inline-block;width:13px;height:13px;border:3px solid #e8a000;border-radius:50%;vertical-align:middle;"></span> core master regulator: significant gene feeding &ge;MINTGT disease genes<br>
 <span style="color:#c1121f;font-weight:bold;">&#8212;&#8212;</span> flux toward risk &nbsp;
 <span style="color:#2c6fbb;font-weight:bold;">&#8212;&#8212;</span> flux away from risk &nbsp;(number on edge = flux)<br>
 <span style="color:#555;font-weight:bold;">- - -</span> dashed = regulator on <b>same chromosome</b> as target (distance in Mb shown; small = possible co-expression/cis)
</div>
<script>
function fluxSave(name, kind){
 var c=document.getElementsByTagName('canvas')[0];
 var t=document.createElement('canvas');t.width=c.width;t.height=c.height;
 var x=t.getContext('2d');x.fillStyle='#fff';x.fillRect(0,0,t.width,t.height);x.drawImage(c,0,0);
 if(kind==='png'){var a=document.createElement('a');a.href=t.toDataURL('image/png');a.download=name;a.click();return;}
 var b64=t.toDataURL('image/jpeg',0.95).split(',')[1],bin=atob(b64);
 var jpg=new Uint8Array(bin.length);for(var i=0;i<bin.length;i++)jpg[i]=bin.charCodeAt(i);
 var w=t.width,h=t.height;
 var enc=function(s){var u=new Uint8Array(s.length);for(var i=0;i<s.length;i++)u[i]=s.charCodeAt(i)&255;return u;};
 var parts=[],off=[],pos=0;
 function push(p){var u=(typeof p==='string')?enc(p):p;parts.push(u);pos+=u.length;}
 function obj(n,head,stream){off[n]=pos;push(n+' 0 obj\\n'+head);if(stream!==undefined){push('\\nstream\\n');push(stream);push('\\nendstream');}push('\\nendobj\\n');}
 push('%PDF-1.3\\n');
 obj(1,'<</Type/Catalog/Pages 2 0 R>>');
 obj(2,'<</Type/Pages/Kids[3 0 R]/Count 1>>');
 obj(3,'<</Type/Page/Parent 2 0 R/MediaBox[0 0 '+w+' '+h+']/Resources<</XObject<</Im0 4 0 R>>>>/Contents 5 0 R>>');
 obj(4,'<</Type/XObject/Subtype/Image/Width '+w+'/Height '+h+'/ColorSpace/DeviceRGB/BitsPerComponent 8/Filter/DCTDecode/Length '+jpg.length+'>>',jpg);
 var ct='q\\n'+w+' 0 0 '+h+' 0 0 cm\\n/Im0 Do\\nQ';
 obj(5,'<</Length '+ct.length+'>>',ct);
 var xp=pos,nn=6,xr='xref\\n0 '+nn+'\\n0000000000 65535 f \\n';
 for(var j=1;j<nn;j++){var s=''+off[j];while(s.length<10)s='0'+s;xr+=s+' 00000 n \\n';}
 push(xr);push('trailer\\n<</Size '+nn+'/Root 1 0 R>>\\nstartxref\\n'+xp+'\\n%%EOF');
 var tot=parts.reduce(function(a,p){return a+p.length;},0),out=new Uint8Array(tot),o=0;
 parts.forEach(function(p){out.set(p,o);o+=p.length;});
 var a=document.createElement('a');a.href=URL.createObjectURL(new Blob([out],{type:'application/pdf'}));a.download=name;a.click();
}
</script>
<div style="position:absolute;top:12px;right:12px;z-index:1000;">
 <button onclick="fluxSave('SHORT_flux_network.pdf','pdf')" style="padding:9px 15px;font-family:Arial;font-size:14px;
  background:#6a51a3;color:#fff;border:none;border-radius:6px;cursor:pointer;box-shadow:0 1px 4px rgba(0,0,0,.3);">&#128196; Save PDF</button>
 <button onclick="fluxSave('SHORT_flux_network.png','png')" style="padding:9px 12px;font-family:Arial;font-size:13px;margin-left:6px;
  background:#888;color:#fff;border:none;border-radius:6px;cursor:pointer;box-shadow:0 1px 4px rgba(0,0,0,.3);">&#128247; PNG</button>
</div>
"""


def build(key, flux_dir, trait, out_dir, S, chrom, mid, min_targets):
    N = pd.read_csv(f"{flux_dir}/flux_nodes_{key}_{trait}.tsv", sep="\t").set_index("gene_id")
    E = pd.read_csv(f"{flux_dir}/flux_edges_{key}_{trait}.tsv", sep="\t")
    if E.empty:
        print(f"[{key}] no flux edges -- skipped"); return
    sig = set(N.index[N["disease_sig"] == True])
    # whole flux map (disease genes + their up-to-2-hop ancestors), largest connected component.
    # A regulator is a hub (highlighted) only if it feeds >= min_targets disease genes:
    #   significant hub -> core master regulator (gold ring);  non-significant hub -> key driver.
    G = nx.from_pandas_edgelist(E, "regulator", "target", edge_attr="flux", create_using=nx.DiGraph())
    G = G.subgraph(max(nx.weakly_connected_components(G), key=len)).copy()
    fmax = max(abs(d["flux"]) for *_, d in G.edges(data=True))
    is_dg = lambda n: n in sig

    net = Network(height="100vh", width="100vw", directed=True, bgcolor="#ffffff",
                  font_color="#222", cdn_resources="in_line")
    net.force_atlas_2based(gravity=-45, central_gravity=0.013, spring_length=95, spring_strength=0.08,
                           damping=0.5, overlap=0.7)
    net.toggle_physics(True)

    for n in G.nodes():
        c = chrom(n)
        ring = is_dg(n) and G.out_degree(n) >= min_targets
        if is_dg(n):
            zc = float(N["z_cis"].get(n, 0.0)); zt = float(N["z_trans"].get(n, 0.0))
            net.add_node(n, label=f"{S(n)}  (chr{c})\nZcis {zc:+.1f} · Ztr {zt:+.1f}",
                         color={"background": C_CORE, "border": C_DG if ring else "#4a2f7a"},
                         size=20 + min(24, abs(zt) * 4), shape="dot", borderWidth=5 if ring else 1.5,
                         title=f"{S(n)} — disease gene · chr{c}\nZ_cis={zc:+.2f}  Z_trans={zt:+.2f}",
                         font={"size": 20, "bold": True, "color": "#222"})
        elif G.out_degree(n) >= min_targets:
            # key driver: non-significant regulator feeding >=min_targets disease genes; size ~ # targets
            od = G.out_degree(n)
            net.add_node(n, label=f"{S(n)}  (chr{c})", color={"background": C_KD, "border": "#0b5e59"},
                         size=13 + min(30, od * 1.4), shape="dot", borderWidth=2,
                         title=f"{S(n)} — key driver · feeds {od} disease genes · chr{c}",
                         font={"size": 14})
        else:
            net.add_node(n, label=f"{S(n)}  (chr{c})", color={"background": C_PRIV, "border": "#6baed6"},
                         size=7, shape="dot", borderWidth=1, title=f"{S(n)} — regulator · chr{c}",
                         font={"size": 10, "color": "#888"})

    for u, v, d in G.edges(data=True):
        f = d["flux"]; col = C_POS if f >= 0 else C_NEG
        cu, cv, mu, mv = chrom(u), chrom(v), mid(u), mid(v)
        same = (cu == cv and cu != "?" and mu is not None and mv is not None)
        if same:                                          # regulator on same chromosome as its target
            dist = abs(mu - mv) / 1e6                      # -> show genomic distance (co-expression clue)
            lab = f"{f:+.2g} · {dist:.2f}Mb"
            note = (f"\nSAME chromosome (chr{cu}) — {dist:.2f} Mb apart"
                    + ("  ⚠ likely co-expression/cis" if dist < 1 else ""))
        else:
            lab = f"{f:+.2g}"
            note = f"\ntrans: regulator chr{cu} → gene chr{cv}"
        net.add_edge(u, v, color=col, width=1 + 6 * (abs(f) / fmax) ** 0.5, label=lab,
                     arrows="to", smooth={"type": "dynamic"}, dashes=same,
                     font={"size": 11, "color": col, "strokeWidth": 4, "strokeColor": "#ffffff", "align": "middle"},
                     title=f"{S(u)} → {S(v)}  flux={f:+.3g}" + note)

    ndg = sum(1 for n in G if is_dg(n))
    label = f"{TISSUE_NAME.get(key, key)} ({key})"
    legend = (LEGEND.replace("LABEL", label).replace("NDG", str(ndg))
              .replace("NE", str(G.number_of_edges())).replace("SHORT", key)
              .replace("MINTGT", str(min_targets)))
    out = f"{out_dir}/fig_flux_lcc_{key}.html"
    net.write_html(out, notebook=False)
    Path(out).write_text(Path(out).read_text().replace("<body>", "<body>\n" + legend, 1))
    print(f"[{key}] nodes {G.number_of_nodes()} ({ndg} disease), edges {G.number_of_edges()} -> {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--flux-dir", required=True,
                    help="dir with flux_nodes_<key>_<trait>.tsv and flux_edges_<key>_<trait>.tsv")
    ap.add_argument("--gencode", required=True, help="GENCODE genes tsv (gene_id, gene_name, chr, start, end)")
    ap.add_argument("--outdir", required=True, help="output dir for the HTML files")
    ap.add_argument("--trait", default="CAD_aragam2022", help="trait label in the flux filenames")
    ap.add_argument("--keys", nargs="+", default=None,
                    help="tissue keys; default: discovered from --flux-dir")
    ap.add_argument("--min-targets", type=int, default=10,
                    help="min disease genes a regulator must feed to be highlighted as a hub "
                         "(key driver if non-significant, core master regulator if significant)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    L = gl.load_label_map(gl.default_annot(args.gencode))
    S = lambda g: L.get(g, L.get(str(g).split(".")[0], str(g).split(".")[0]))
    gc = pd.read_csv(args.gencode, sep="\t")
    gc["base"] = gc["gene_id"].str.split(".").str[0]
    CHR = gc.set_index("base")["chr"].astype(str).to_dict()
    MID = (((gc["start"] + gc["end"]) / 2).groupby(gc["base"]).first()).to_dict()
    chrom = lambda g: CHR.get(str(g).split(".")[0], "?")
    mid = lambda g: MID.get(str(g).split(".")[0])

    keys = args.keys or sorted(
        os.path.basename(f)[len("flux_nodes_"):].rsplit(f"_{args.trait}.tsv", 1)[0]
        for f in glob.glob(f"{args.flux_dir}/flux_nodes_*_{args.trait}.tsv"))
    if not keys:
        raise SystemExit(f"no flux_nodes_*_{args.trait}.tsv found in {args.flux_dir}")
    for key in keys:
        build(key, args.flux_dir, args.trait, args.outdir, S, chrom, mid, args.min_targets)


if __name__ == "__main__":
    main()
