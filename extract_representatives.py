"""Extract one representative protein per Panaroo family for eggNOG-mapper.
Usage: python3 -I extract_representatives.py <panaroo_dir> <out.faa>
"""
import re, sys
from pathlib import Path
import pandas as pd

panaroo, out = Path(sys.argv[1]), Path(sys.argv[2])
gene = pd.read_csv(panaroo / "gene_data.csv", dtype=str).fillna("")
seq = dict(zip(gene["clustering_id"], gene["prot_sequence"]))
graph = (panaroo / "final_graph.gml").read_text(errors="ignore")
n = 0
with open(out, "w") as fh:
    for block in re.findall(r'node \[(.*?)\n  \]', graph, re.S):
        fam = re.search(r'name "([^"]+)"', block)
        cen = re.search(r'longCentroidID "([^"]+)"', block)
        if not fam or not cen:
            continue
        rep = next((c for c in cen.group(1).split(";") if c in seq), None)
        if rep is None:
            continue
        s = seq[rep]
        fh.write(f">{fam.group(1)}|{rep}\n")
        fh.write("\n".join(s[i:i + 80] for i in range(0, len(s), 80)) + "\n")
        n += 1
print(f"{n} representative sequences written to {out}")
