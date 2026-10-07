"""Re-analysis of the Shewanella benthica Panaroo/eggNOG pangenome tables.
Usage: python3 analysis.py <input.xlsx> <outdir>
"""
import sys, re, json, itertools
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, chi2_contingency
from scipy.optimize import curve_fit
def multipletests(p, method="fdr_bh"):
    p = np.asarray(p, float); n = len(p); o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]; out = np.empty(n); out[o] = np.minimum(q, 1)
    return None, out

xlsx, out = sys.argv[1], sys.argv[2]
X = pd.read_excel(xlsx, sheet_name=None)
STRAINS = ["BOEU19-1", "DB21MT-2", "DSM8812", "KT99"]
HADAL = {"DB21MT-2", "KT99"}
CATS = ["Core", "Accessory", "Strain-specific"]
R = {}

m = X["Master_functional"].copy()
for c in ["Strain", "Preferred_name", "Description", "COG_category", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "PFAMs"]:
    m[c] = m[c].fillna("-").astype(str).str.strip()
m["txt"] = (m.Panaroo_family + " " + m.Preferred_name + " " + m.Description + " " + m.PFAMs).str.lower()

# ---------------- annotation coverage (corrected) ----------------
m["has_eggnog"] = ~((m.Description == "-") & (m.COG_category == "-") & (m.PFAMs == "-") & (m.Preferred_name == "-"))
m["has_desc"] = m.Description != "-"
m["has_cog"] = m.COG_category.str.fullmatch(r"[A-Z]+")
m["has_cog_informative"] = m.has_cog & ~m.COG_category.isin(["S", "R"])
m["has_kegg"] = m.KEGG_ko.str.startswith("ko:")
cov = []
for c in CATS:
    s = m[m.Genome_category == c]
    n = len(s)
    cov.append(dict(Genome_category=c, Families=n,
                    eggNOG_hit=int(s.has_eggnog.sum()), Description=int(s.has_desc.sum()),
                    COG_assigned=int(s.has_cog.sum()), COG_informative=int(s.has_cog_informative.sum()),
                    KEGG_KO=int(s.has_kegg.sum())))
cov = pd.DataFrame(cov)
tot = cov.sum(numeric_only=True); tot["Genome_category"] = "Total"
cov = pd.concat([cov, pd.DataFrame([tot])], ignore_index=True)
for k in ["eggNOG_hit", "Description", "COG_assigned", "COG_informative", "KEGG_KO"]:
    cov[k + "_pct"] = (100 * cov[k] / cov.Families).round(1)
R["coverage"] = cov.to_dict("records")

# ---------------- curated functional classes ----------------
TNP_P = r"dde_tnp|transposase|hth_tnp|\btnp|y1_tnp|y2_tnp|zf-is66|dedd_tnp|dimer_tnp|orfb_is605|orfb_zn_ribbon|isxo2|insertion element|cog1943|cog2801|cog3328|cog3436|cog3464|cog4584|cog2963|cog3039|cog3316|cog3547|cog5421|cog5433|cog4644|\bis[0-9]{1,4}\b|\bisl3\b|\btnpa|\btnpb|\binsa|\binsb|\bistb|mule"
PHAGE_P = r"\bphage\b(?! shock)|prophage|capsid|phage terminase|terminase (large|small)|gpp-like|portal protein|baseplate|tail fiber|tail tube|tail sheath|\bholin\b|integrase|site-specific recombinase|resolvase(?!, rnase h)|excisionase|phage_int|\bxis\b"
PLASMID_P = r"conjugal|conjugative|relaxase|f plasmid|plasmid transfer|traf_2|t4ss|type iv secretion system|mobilization protein|plasmid partition"
CRISPR_P = r"crispr-associated|clustered regularly interspaced|\bcas_|\bcas[0-9]+[a-z]?\b|\bcsy[0-9]\b|\bcse[0-9]\b"
RM_P = r"restriction|\bhsd[rms]|methylase_s|hsdm_n|hsdr|mrr_cat|mvai|n6_mtase|dna methylase|dna adenine methylase|cytosine-specific methyltransferase|dcm\b"
TA_P = r"toxin-antitoxin|antitoxin|(?<!tcdb)_toxin\b|-like toxin|abiei|abieii|abortive infection|\brele\b|\bhig[ab]\b|\bhic[ab]\b|\bmaz[ef]\b|\bvap[bc]\b|\byafq\b|\byafo\b|\bdinj\b|\byoeb\b|\bphd\b|mqsa"

EXCL = {"xerc", "xerd", "rnc", "vgrg", "pare", "cmr", "dnaj"}
def excluded(r):
    return r.Preferred_name.lower() in EXCL
def classify(r):
    t = r.txt
    if excluded(r): t = ""
    if re.search(TNP_P, t):
        return "Mobilome: IS elements/transposases"
    if re.search(PHAGE_P, t) or re.search(PLASMID_P, t):
        return "Mobilome: prophages, integrases and plasmids"
    if re.search(CRISPR_P, t) or re.search(RM_P, t) or re.search(TA_P, t):
        return "Defence systems (CRISPR-Cas, R-M, TA)"
    cog = r.COG_category
    if not re.fullmatch(r"[A-Z]+", cog) or cog[0] in "SR":
        return "Unknown / poorly characterized"
    c = cog[0]
    if c == "M":
        return "Cell envelope biogenesis (M)"
    if c in "TK":
        return "Signal transduction and transcription (T, K)"
    if c in "PU":
        return "Inorganic ion transport and secretion (P, U)"
    if c in "CEFGHIQ":
        return "Metabolism (C, E, F, G, H, I, Q)"
    return "Other characterized functions"

def subclass(r):
    t = r.txt
    if excluded(r): return ""
    if re.search(CRISPR_P, t): return "CRISPR-Cas"
    if re.search(RM_P, t): return "Restriction-modification"
    if re.search(TA_P, t): return "Toxin-antitoxin / abortive infection"
    if re.search(TNP_P, t): return "Transposase / IS element"
    if re.search(r"integrase|recombinase|resolvase|excisionase|phage_int|\bxis\b", t): return "Integrase / recombinase"
    if re.search(PHAGE_P, t): return "Phage structural / lysis"
    if re.search(PLASMID_P, t): return "Plasmid / conjugation"
    return ""

m["Functional_class"] = m.apply(classify, axis=1)
m["Mobilome_defence_subclass"] = m.apply(subclass, axis=1)
CLASS_ORDER = ["Mobilome: IS elements/transposases", "Mobilome: prophages, integrases and plasmids",
               "Defence systems (CRISPR-Cas, R-M, TA)", "Cell envelope biogenesis (M)",
               "Signal transduction and transcription (T, K)", "Inorganic ion transport and secretion (P, U)",
               "Metabolism (C, E, F, G, H, I, Q)", "Other characterized functions", "Unknown / poorly characterized"]

# IS family assignment from Pfam domain / description
IS_RULES = [("IS66", r"is66|zf-is66|dde_tnp_is66"), ("IS4/IS5 (DDE_Tnp_1)", r"dde_tnp_1\b|dde_tnp_1,|dde_tnp_1$|dimer_tnp_tn5|is4|is5\b|tnp_dna_bind"),
            ("IS3/IS911", r"is3|is911|hth_tnp_1|rve_3|hth_21|cog2801|cog2963"), ("IS21", r"is21|istb|cog4584"),
            ("IS110", r"is110|dedd_tnp_is110|transposase_20"), ("IS200/IS605", r"is200|is605|y1_tnp|orfb|cog1943"),
            ("IS91", r"is91|y2_tnp"), ("IS1595/ISXO2", r"is1595|isxo2"), ("IS256/MULE", r"is256|transposase_mut|mule"),
            ("IS30/IS630", r"is30|is630|hth_tnp_tc|dde_3")]
def isfam(t):
    for k, p in IS_RULES:
        if re.search(p, t): return k
    return "Unassigned transposase"
m["IS_family"] = np.where(m.Mobilome_defence_subclass == "Transposase / IS element", m.txt.map(isfam), "")

# ---------------- class composition ----------------
comp = pd.crosstab(m.Functional_class, m.Genome_category).reindex(CLASS_ORDER).reindex(columns=CATS).fillna(0).astype(int)
ss = m[m.Genome_category == "Strain-specific"]
comp_strain = pd.crosstab(ss.Functional_class, ss.Strain).reindex(CLASS_ORDER).reindex(columns=STRAINS).fillna(0).astype(int)
R["class_by_compartment"] = comp.to_dict()
R["class_by_strain"] = comp_strain.to_dict()

sub_strain = pd.crosstab(ss.Mobilome_defence_subclass, ss.Strain).drop(index="", errors="ignore").reindex(columns=STRAINS).fillna(0).astype(int)
sub_comp = pd.crosstab(m.Mobilome_defence_subclass, m.Genome_category).drop(index="", errors="ignore").reindex(columns=CATS).fillna(0).astype(int)
is_strain = pd.crosstab(ss[ss.IS_family != ""].IS_family, ss[ss.IS_family != ""].Strain).reindex(columns=STRAINS).fillna(0).astype(int)
R["subclass_by_strain"] = sub_strain.to_dict(); R["subclass_by_compartment"] = sub_comp.to_dict(); R["IS_by_strain"] = is_strain.to_dict()

# mobilome tests
mob = m.Functional_class.str.startswith("Mobilome")
ism = m.Functional_class == "Mobilome: IS elements/transposases"
tests = []
def fisher_row(name, a, b, c, d, desc):
    o, p = fisher_exact([[a, b], [c, d]])
    tests.append(dict(Test=name, Description=desc, a=a, b=b, c=c, d=d, odds_ratio=o, p_value=p))
for lab, mask in [("Mobilome (all)", mob), ("IS elements/transposases", ism)]:
    sS = (m.Genome_category == "Strain-specific"); sC = (m.Genome_category == "Core"); sA = (m.Genome_category == "Accessory")
    fisher_row(f"{lab}: strain-specific vs core", int((mask & sS).sum()), int((~mask & sS).sum()), int((mask & sC).sum()), int((~mask & sC).sum()), "[[SS with, SS without],[Core with, Core without]]")
    fisher_row(f"{lab}: accessory vs core", int((mask & sA).sum()), int((~mask & sA).sum()), int((mask & sC).sum()), int((~mask & sC).sum()), "[[Acc with, Acc without],[Core with, Core without]]")
    h = sS & m.Strain.isin(HADAL); nh = sS & ~m.Strain.isin(HADAL)
    fisher_row(f"{lab}: hadal vs non-hadal strain-specific", int((mask & h).sum()), int((~mask & h).sum()), int((mask & nh).sum()), int((~mask & nh).sum()), "[[hadal SS with, without],[non-hadal SS with, without]]")
defm = m.Functional_class.str.startswith("Defence")
sS = (m.Genome_category == "Strain-specific"); h = sS & m.Strain.isin(HADAL); nh = sS & ~m.Strain.isin(HADAL)
fisher_row("Defence: hadal vs non-hadal strain-specific", int((defm & h).sum()), int((~defm & h).sum()), int((defm & nh).sum()), int((~defm & nh).sum()), "[[hadal SS with, without],[non-hadal SS with, without]]")
unk = m.Functional_class.str.startswith("Unknown")
fisher_row("Unknown: strain-specific vs core", int((unk & sS).sum()), int((~unk & sS).sum()), int((unk & (m.Genome_category == 'Core')).sum()), int((~unk & (m.Genome_category == 'Core')).sum()), "[[SS with, without],[Core with, without]]")
tests = pd.DataFrame(tests)
tests["p_BH"] = multipletests(tests.p_value, method="fdr_bh")[1]
R["fisher_tests"] = tests.to_dict("records")
chi = chi2_contingency(comp.values)
R["chi2_class_x_compartment"] = dict(chi2=chi[0], p=chi[1], dof=chi[2])
chi_s = chi2_contingency(comp_strain.values)
R["chi2_class_x_strain"] = dict(chi2=chi_s[0], p=chi_s[1], dof=chi_s[2])

# ---------------- COG enrichment (families counted once per letter) ----------------
LET = list("ABCDEFGHIJKLMNOPQSTUVWZ")
cogm = m[m.has_cog]
rows = []
for L in LET:
    rec = {"COG": L}
    for c in CATS:
        s = cogm[cogm.Genome_category == c]
        rec[f"{c}_n"] = int(s.COG_category.str.contains(L).sum()); rec[f"{c}_N"] = len(s)
        rec[f"{c}_pct"] = 100 * rec[f"{c}_n"] / len(s)
    for c in ["Accessory", "Strain-specific"]:
        a, b = rec[f"{c}_n"], rec[f"{c}_N"] - rec[f"{c}_n"]
        cc, d = rec["Core_n"], rec["Core_N"] - rec["Core_n"]
        o, p = fisher_exact([[a, b], [cc, d]])
        rec[f"{c}_vs_Core_OR"] = o; rec[f"{c}_vs_Core_p"] = p
        rec[f"{c}_vs_Core_log2FC"] = np.log2(((a + .5) / (rec[f'{c}_N'] + 1)) / ((cc + .5) / (rec['Core_N'] + 1)))
    rows.append(rec)
cog = pd.DataFrame(rows)
cog = cog[(cog[[f"{c}_n" for c in CATS]].sum(axis=1)) > 0].reset_index(drop=True)
for c in ["Accessory", "Strain-specific"]:
    cog[f"{c}_vs_Core_q"] = multipletests(cog[f"{c}_vs_Core_p"], method="fdr_bh")[1]
R["cog"] = cog.to_dict("records")

# ---------------- intersections -> exact subset rarefaction ----------------
inter = X["Intersections"].copy()
inter["set"] = inter.Combination.str.split(r" \+ ").map(frozenset)
fam_per = {s: int(inter[inter.set.map(lambda x: s in x)].Families.sum()) for s in STRAINS}
R["families_per_strain"] = fam_per
pair = {}
for a, b in itertools.combinations(STRAINS, 2):
    sh = int(inter[inter.set.map(lambda x: a in x and b in x)].Families.sum())
    un = int(inter[inter.set.map(lambda x: a in x or b in x)].Families.sum())
    excl = int(inter[inter.set.map(lambda x: x == frozenset([a, b]))].Families.sum())
    pair[f"{a}|{b}"] = dict(shared=sh, union=un, jaccard=sh / un, exclusive_pair=excl)
R["pairwise"] = pair
sub = []
for k in range(1, 5):
    for S in itertools.combinations(STRAINS, k):
        S = set(S)
        pan = int(inter[inter.set.map(lambda x: len(x & S) > 0)].Families.sum())
        core = int(inter[inter.set.map(lambda x: S <= x)].Families.sum())
        sub.append(dict(n_genomes=k, strains=" + ".join(sorted(S)), pangenome=pan, core=core))
sub = pd.DataFrame(sub)
heaps = lambda N, K, g: K * N ** g
popt, pcov = curve_fit(heaps, sub.n_genomes, sub.pangenome, p0=[3600, .3])
perr = np.sqrt(np.diag(pcov))
pred = heaps(sub.n_genomes, *popt); ss_res = ((sub.pangenome - pred) ** 2).sum(); ss_tot = ((sub.pangenome - sub.pangenome.mean()) ** 2).sum()
R["heaps_panaroo"] = dict(K=popt[0], gamma=popt[1], K_se=perr[0], gamma_se=perr[1], R2=1 - ss_res / ss_tot)
means = sub.groupby("n_genomes")[["pangenome", "core"]].agg(["mean", "std", "min", "max"])
R["rarefaction_means"] = {int(k): {f"{a}_{b}": float(v) for (a, b), v in row.items()} for k, row in means.iterrows()}
# new genes added by the k-th genome (mean over all orderings)
newg = []
for perm in itertools.permutations(STRAINS):
    seen = set()
    for i in range(4):
        before = int(inter[inter.set.map(lambda x: len(x & seen) > 0)].Families.sum()) if seen else 0
        seen = seen | {perm[i]}
        after = int(inter[inter.set.map(lambda x: len(x & seen) > 0)].Families.sum())
        newg.append(dict(k=i + 1, new=after - before))
newg = pd.DataFrame(newg).groupby("k").new.mean()
R["new_genes_mean"] = newg.to_dict()
pw = lambda N, k, a: k * N ** (-a)
po, _ = curve_fit(pw, newg.index[1:].values.astype(float), newg.values[1:], p0=[1500, 1])
R["new_gene_powerlaw"] = dict(kappa=po[0], alpha=po[1])

# ---------------- adaptive trait inventory ----------------
# Trait inventory: gene symbols matched on names (Panaroo family + eggNOG preferred name, "_n" suffixes stripped),
# phrases matched on description + Pfam. Each tuple: (module, trait, symbol_regex, phrase_regex)
m["names"] = (m.Panaroo_family.str.replace(r"_[0-9]+", "", regex=True).str.replace("~~~", " ") + " " + m.Preferred_name).str.lower()
m["descp"] = (m.Description + " " + m.PFAMs).str.lower()
TRAITS = [
 ("Membrane homeoviscosity", "PUFA (EPA) synthase PfaA/B/C/D", r"\bpfa[a-e]\b", r"polyunsaturated fatty acid synth"),
 ("Membrane homeoviscosity", "Fatty acid desaturase, FabF, FabB", r"\bdes[a-c]\b|\bfab[fb]\b", r"fatty acid desaturase"),
 ("Respiration & bioenergetics", "Na+-NQR (NqrA-F, NqrM)", r"\bnqr[a-fm]\b", r"nqr complex|nqrm"),
 ("Respiration & bioenergetics", "Rnf complex", r"\brnf[a-g]\b", r"^$x"),
 ("Respiration & bioenergetics", "Cytochrome bd oxidase (CydAB)", r"\bcyd[ab]\b", r"cytochrome bd"),
 ("Respiration & bioenergetics", "TMAO reductase system (TorA/C/D/Y, TorRST)", r"\btor[acdrsty]\b", r"trimethylamine"),
 ("Respiration & bioenergetics", "Nitrate/nitrite reduction (Nap, Nrf, Nir)", r"\bnap[a-h]\b|\bnrf[a-h]\b|\bnir[bdk]\b", r"nitrate reductase|nitrite reductase"),
 ("Respiration & bioenergetics", "Multiheme c-type cytochromes (Mtr/Omc)", r"\bmtr[abc]\b|\bomc[a-z]\b", r"decaheme|multiheme"),
 ("Piezolytes & polyamines", "Polyamine biosynthesis (SpeA/B/D/E)", r"\bspe[abde]\b", r"spermidine synthase|spermine"),
 ("Piezolytes & polyamines", "Polyamine ABC transport (PotABCD/FGHI)", r"\bpot[a-i]\b", r"^$x"),
 ("Piezolytes & polyamines", "Glycine betaine / choline / ectoine", r"\bbet[abt]\b|\bect[abcd]\b|\bopu[a-d]\b|\bpro[vwx]\b|\bbcct\b", r"glycine betaine|ectoine|bcct transporter|bcct_"),
 ("Piezolytes & polyamines", "Mechanosensitive channels (MscS/MscL)", r"\bmsc[sl]\b", r"mechanosensitive"),
 ("Stress & genome maintenance", "Chaperones/proteases (DnaK/J, GroESL, Clp, HtpG, Lon)", r"\bdna[kj]\b|\bgro[els]+\b|\bclp[abpxs]\b|\bhtp[gx]\b|\bhsl[uv]\b|\blon\b", r"^$x"),
 ("Stress & genome maintenance", "Cold-shock proteins (Csp)", r"\bcsp[a-z]\b", r"cold.shock"),
 ("Stress & genome maintenance", "Recombination and SOS repair (Rec, Uvr, UmuDC, LexA)", r"\brec[abcdfnoqrjg]\b|\buvr[abcd]\b|\bumu[cd]\b|\blexa\b", r"^$x"),
 ("Surface & motility", "Flagellar apparatus", r"\bfl[gih][a-z]\b|\bmot[ab]\b", r"flagell"),
 ("Surface & motility", "Type IV pili", r"\bpil[a-z]\b", r"type iv pil|pilin"),
 ("Surface & motility", "Capsule/LPS/O-antigen biosynthesis", r"\bwz[a-z]\b|\bwec[a-g]\b|\brf[abce]\b|\brml[a-d]\b|\bkps[a-z]\b|\bwb[a-z]{2}\b|\bwf[a-z]{2}\b", r"lipopolysacch|exopolysacch|capsular polysacch|o-antigen"),
 ("Nutrient acquisition", "Chitinases and glycoside hydrolases", r"\bchi[abc]\b", r"chitin|glyco_hydro|glycosyl hydrolase|glycoside hydrolase"),
 ("Nutrient acquisition", "TonB-dependent iron/siderophore uptake", r"\btonb\b|\bfhu[a-e]\b|\bfep[a-g]\b|\bfec[a-e]\b", r"tonb|siderophore|ferric|hemin|heme receptor"),
 ("Nutrient acquisition", "Heavy-metal efflux/resistance (Czc, Cus, Cop, Mer, Ars)", r"\bczc[abcd]\b|\bcus[abcf]\b|\bcop[ab]\b|\bmer[aprt]\b|\bars[abcr]\b", r"cobalt-zinc-cadmium|copper resistance|mercur|arsenate|arsenical"),
 ("Mobilome & defence", "Transposases / IS elements", r"^$x", None),
 ("Mobilome & defence", "Prophages, integrases and plasmid functions", r"^$x", None),
 ("Mobilome & defence", "CRISPR-Cas (type I-F, Csy)", r"^$x", None),
 ("Mobilome & defence", "Restriction-modification", r"^$x", None),
 ("Mobilome & defence", "Toxin-antitoxin / abortive infection", r"^$x", None),
]
SUBMAP = {"Transposases / IS elements": ["Transposase / IS element"],
          "Prophages, integrases and plasmid functions": ["Integrase / recombinase", "Phage structural / lysis", "Plasmid / conjugation"],
          "CRISPR-Cas (type I-F, Csy)": ["CRISPR-Cas"], "Restriction-modification": ["Restriction-modification"],
          "Toxin-antitoxin / abortive infection": ["Toxin-antitoxin / abortive infection"]}
tr = []
m["Trait"] = ""
for grp, name, sym, phr in TRAITS:
    if phr is None:
        mask = m.Mobilome_defence_subclass.isin(SUBMAP[name])
    else:
        mask = m.names.str.contains(sym, regex=True) | m.descp.str.contains(phr, regex=True)
    s = m[mask]
    m.loc[mask & (m.Trait == ""), "Trait"] = name
    rec = dict(Module=grp, Trait=name, Total=len(s))
    for c in CATS: rec[c] = int((s.Genome_category == c).sum())
    for st in STRAINS: rec[f"SS_{st}"] = int(((s.Genome_category == "Strain-specific") & (s.Strain == st)).sum())
    rec["Example_genes"] = ", ".join(sorted(set([x for x in s.Preferred_name if x != "-"]))[:16])
    tr.append(rec)
traits = pd.DataFrame(tr)
R["traits"] = traits.to_dict("records")

# ---------------- save ----------------
import os; os.makedirs(out, exist_ok=True)
m.drop(columns=["txt","names","descp"]).to_pickle(f"{out}/master.pkl")
comp.to_csv(f"{out}/class_by_compartment.csv"); comp_strain.to_csv(f"{out}/class_by_strain.csv")
sub_strain.to_csv(f"{out}/subclass_by_strain.csv"); sub_comp.to_csv(f"{out}/subclass_by_compartment.csv"); is_strain.to_csv(f"{out}/IS_by_strain.csv")
tests.to_csv(f"{out}/fisher_tests.csv", index=False); cog.to_csv(f"{out}/cog_enrichment.csv", index=False)
sub.to_csv(f"{out}/subset_rarefaction.csv", index=False); traits.to_csv(f"{out}/traits.csv", index=False); cov.to_csv(f"{out}/coverage.csv", index=False)
json.dump(R, open(f"{out}/results.json", "w"), indent=1, default=float)
print(json.dumps({k: R[k] for k in ["families_per_strain", "heaps_panaroo", "new_genes_mean", "new_gene_powerlaw", "chi2_class_x_compartment", "chi2_class_x_strain"]}, indent=1, default=float))
print(cov.to_string()); print(comp); print(comp_strain); print(sub_strain); print(sub_comp); print(is_strain)
print(tests.to_string()); print(pd.DataFrame(pair).T)
print(cog[["COG", "Core_pct", "Accessory_pct", "Strain-specific_pct", "Strain-specific_vs_Core_log2FC", "Strain-specific_vs_Core_q", "Accessory_vs_Core_log2FC", "Accessory_vs_Core_q"]].round(4).to_string())
print(traits.drop(columns="Example_genes").to_string())
print(sub.to_string())
