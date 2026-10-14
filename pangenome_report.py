"""Build the family-level workbook from Panaroo and eggNOG-mapper output.

Merges the Panaroo presence/absence table with the eggNOG-mapper annotation of
one representative protein per family (see extract_representatives.py) and
writes pangenome_results.xlsx, the input of analysis.py and figures.py.

Usage:
    python3 pangenome_report.py <gene_presence_absence.csv> <eggnog.emapper.annotations> <out.xlsx>
"""
import sys
from itertools import combinations

import numpy as np
import pandas as pd

# Panaroo column name -> short strain name (order used throughout the workbook)
GENOMES = {
    "Shewanella_benthica_BOEU19-1_prokka": "BOEU19-1",
    "Shewanella_benthica_DB21MT-2_prokka": "DB21MT-2",
    "Shewanella_benthica_DSM8812_prokka": "DSM8812",
    "Shewanella_benthica_KT99_prokka": "KT99",
}
CATS = ["Core", "Accessory", "Strain-specific"]
COG_ORDER = list("ABCDEFGHIKLMNOPQRSTUV")
EGG_COLS = [
    "query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",
    "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko",
    "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE",
    "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs",
]


def valid(x):
    return str(x).strip() not in {"", "-", "nan", "None"}


def load(pa_file, egg_file):
    pa = pd.read_csv(pa_file, dtype=str).fillna("")
    missing = [g for g in GENOMES if g not in pa.columns]
    if missing:
        sys.exit("Missing genome columns in presence/absence table: " + ", ".join(missing))

    present = pd.DataFrame({s: pa[g].map(valid) for g, s in GENOMES.items()})
    pa["n_strains"] = present.sum(axis=1).astype(int)
    pa["Genome_category"] = np.select(
        [pa.n_strains.eq(len(GENOMES)), pa.n_strains.eq(1)], ["Core", "Strain-specific"], default="Accessory")
    pa["Strain"] = np.where(pa.n_strains.eq(1), present.idxmax(axis=1), "")

    egg = pd.read_csv(egg_file, sep="\t", comment="#", header=None, names=EGG_COLS, dtype=str).fillna("")
    # FASTA headers are ">family|representative_id"
    egg["Panaroo_family"] = egg["query"].str.split("|").str[0].str.strip()
    # Safeguard: keep a single record per family, preferring the most informative one
    egg["_score"] = sum(egg[c].map(valid).astype(int) for c in ["Description", "Preferred_name", "COG_category", "KEGG_ko"])
    egg = egg.sort_values(["Panaroo_family", "_score"], ascending=[True, False]).drop_duplicates("Panaroo_family")

    master = (pa[["Gene", "n_strains", "Genome_category", "Strain"]]
              .rename(columns={"Gene": "Panaroo_family"})
              .merge(egg.drop(columns="_score"), on="Panaroo_family", how="left"))
    return pa, present, master


def pangenome_tables(pa, present):
    n = len(pa)
    comp = pa.Genome_category.value_counts().reindex(CATS, fill_value=0).rename_axis("Genome_category").reset_index(name="Families")
    comp["Percentage"] = comp.Families / n * 100
    freq = pa.n_strains.value_counts().sort_index().rename_axis("Number_of_strains").reset_index(name="Families")
    freq["Percentage"] = freq.Families / n * 100
    ss = (pa[pa.Genome_category == "Strain-specific"].Strain.value_counts()
          .reindex(list(GENOMES.values()), fill_value=0).rename_axis("Strain").reset_index(name="Strain_specific_families"))
    strains = list(GENOMES.values())
    rows = []
    for r in range(1, len(strains) + 1):
        for combo in combinations(strains, r):
            mask = present[list(combo)].all(axis=1) & ~present[[s for s in strains if s not in combo]].any(axis=1)
            rows.append(dict(Combination=" + ".join(combo), n_strains=r, Families=int(mask.sum())))
    inter = pd.DataFrame(rows).sort_values(["n_strains", "Families"], ascending=[False, False])
    return comp, freq, ss, inter


def annotation_tables(master):
    cov = []
    for c in CATS:
        s = master[master.Genome_category == c]
        rec = dict(Genome_category=c, Families=len(s))
        for lab, col in [("With_description", "Description"), ("With_COG", "COG_category"), ("With_KEGG", "KEGG_ko")]:
            rec[lab] = int(s[col].map(valid).sum())
        cov.append(rec)
    cov = pd.DataFrame(cov)
    for lab, col in [("Description_pct", "With_description"), ("COG_pct", "With_COG"), ("KEGG_pct", "With_KEGG")]:
        cov[lab] = (cov[col] / cov.Families.replace(0, np.nan) * 100).fillna(0).round(2)

    # COG letters (multi-letter assignments counted once per letter)
    x = master.loc[master.COG_category.map(valid), ["Panaroo_family", "Genome_category", "COG_category"]].copy()
    x["COG_category"] = x.COG_category.str.findall(r"[A-Z]")
    x = x.explode("COG_category")
    cog = x.groupby(["Genome_category", "COG_category"]).size().reset_index(name="Families")
    cog["Percentage"] = cog.Families / cog.groupby("Genome_category").Families.transform("sum") * 100
    heat = (cog.pivot(index="COG_category", columns="Genome_category", values="Percentage")
            .reindex(index=COG_ORDER, columns=CATS).fillna(0).round(3))
    return cov, cog.round(3), heat


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    pa_file, egg_file, out = sys.argv[1:]
    pa, present, master = load(pa_file, egg_file)
    comp, freq, ss, inter = pangenome_tables(pa, present)
    cov, cog, heat = annotation_tables(master)

    qc = pd.DataFrame({
        "Metric": ["Total Panaroo gene families", "Number of strains", "Families with eggNOG description",
                   "Families with COG category", "Families with KEGG KO"],
        "Value": [len(pa), len(GENOMES), int(master.Description.map(valid).sum()),
                  int(master.COG_category.map(valid).sum()), int(master.KEGG_ko.map(valid).sum())],
    })
    cols = ["Panaroo_family", "Genome_category", "Strain", "n_strains", "Preferred_name", "Description",
            "COG_category", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "PFAMs"]
    master_out = master[cols]
    top_ss = (master[master.Genome_category == "Strain-specific"]
              [["Panaroo_family", "Strain", "Preferred_name", "Description", "COG_category", "KEGG_ko", "KEGG_Pathway", "PFAMs"]]
              .sort_values(["Strain", "Preferred_name"], na_position="last"))

    with pd.ExcelWriter(out, engine="openpyxl") as w:
        for df, name in [(qc, "QC_summary"), (comp, "Pangenome_compartments"), (freq, "Family_frequency"),
                         (ss, "Strain_specific"), (inter, "Intersections"), (cov, "Functional_summary"),
                         (cog, "COG_summary"), (heat.reset_index(), "COG_heatmap"),
                         (top_ss, "Top_strain_specific"), (master_out, "Master_functional")]:
            df.to_excel(w, sheet_name=name, index=False)

    print(f"Gene families: {len(pa):,}")
    for _, r in comp.iterrows():
        print(f"  {r.Genome_category:16s} {r.Families:6,d} ({r.Percentage:.2f}%)")
    print("Strain-specific:", dict(zip(ss.Strain, ss.Strain_specific_families)))
    print(cov[["Genome_category", "Families", "Description_pct", "KEGG_pct"]].to_string(index=False))
    print(f"Written: {out}")


if __name__ == "__main__":
    main()
