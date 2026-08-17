# -*- coding: utf-8 -*-
"""Summarize PheWAS results for the 40 coffee instrument SNPs -> Supplementary Table S15."""
import io, sys, re
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
import pandas as pd

RAW = r"C:\Users\曹泽众\Documents\kimi\workspace\EJPC_修订产出\stats\phewas_coffee_snps_raw.csv"
df = pd.read_csv(RAW)
print("rows:", len(df), "| unique traits:", df['trait'].nunique(), "| SNPs with >=1 hit:", df['rsid'].nunique())

CATEGORIES = {
    "Adiposity / body composition": r"body mass index|BMI|weight|waist|hip|fat|lean mass|impedance|basal metabolic|obesity|adipos|trunk|arm |leg |whole body|appendicular|body fat",
    "Coffee / caffeine related": r"coffee|caffeine",
    "Glycaemic / diabetes": r"diabet|glucose|glyc|insulin|HbA1c",
    "Lipids": r"cholesterol|triglycerid|lipid|HDL|LDL|lipoprotein|ApoA|ApoB",
    "Blood pressure": r"blood pressure|systolic|diastolic|pulse pressure|hypertension",
    "Smoking / alcohol": r"smok|cigarett|tobacco|alcohol|drink",
    "Blood cell traits": r"platelet|white blood|red blood|red cell|reticulocyte|lymphocyte|monocyte|neutrophil|eosinophil|basophil|haemoglobin|hemoglobin|haematocrit|hematocrit|mean corpuscular|blood cell|myeloid",
    "Kidney / liver function": r"glomerular|eGFR|creatinine|urea|cystatin|albumin|alkaline phosphatase|alanine|aspartate|gamma|bilirubin|liver|kidney",
    "Height / growth": r"\bheight\b|sitting height|standing",
    "Sleep / chronotype": r"sleep|insomnia|chronotype|morning person",
    "Gene expression (eQTL)": r"^ENSG",
}

def categorize(trait):
    for cat, pat in CATEGORIES.items():
        if re.search(pat, str(trait), re.I):
            return cat
    return "Other"

df["category"] = df["trait"].map(categorize)

# ---- category-level summary ----
cat_sum = (df.groupby("category")
             .agg(associations=("rsid", "size"),
                  snps=("rsid", "nunique"),
                  traits=("trait", "nunique"))
             .sort_values("associations", ascending=False))
print("\n=== Category summary ===")
print(cat_sum.to_string())

# ---- per-SNP profile ----
snp_hits = df.groupby("rsid").agg(hits=("trait", "size"), traits=("trait", "nunique"))
snp_hits = snp_hits.sort_values("hits", ascending=False)
print("\n=== Top 15 SNPs by PheWAS hits ===")
print(snp_hits.head(15).to_string())

# ---- key SNPs of interest ----
for s in ["rs1421085", "rs13387939", "rs2107308", "rs6062682"]:
    sub = df[df.rsid == s]
    cats = sub["category"].value_counts()
    print(f"\n--- {s}: {len(sub)} hits ---")
    print(cats.head(8).to_string())

# ---- adiposity associations per SNP (direction) ----
ad = df[df.category == "Adiposity / body composition"]
print("\nAdiposity hits:", len(ad), "across", ad['rsid'].nunique(), "SNPs")

# ---- save supplementary table S15 (compact: SNP x category counts + top traits) ----
piv = df.pivot_table(index="rsid", columns="category", values="trait", aggfunc="count", fill_value=0)
piv["total_hits"] = piv.sum(axis=1)
piv = piv.sort_values("total_hits", ascending=False)
piv.to_csv(r"C:\Users\曹泽众\Documents\kimi\workspace\EJPC_修订产出\stats\phewas_snp_category_matrix.csv")
print("\nsaved phewas_snp_category_matrix.csv")
print(piv.head(8).to_string())

# ---- exports for Supplementary Table S15 ----
cat_out = cat_sum.reset_index()
cat_out.columns = ["category", "associations", "snps", "traits"]
cat_out.to_csv(r"C:\Users\曹泽众\Documents\kimi\workspace\EJPC_修订产出\stats\phewas_category_summary.csv",
               index=False)

ref40 = set(pd.read_csv(r"D:\Desktop\EJPC_投稿_咖啡与房颤\03_Supplementary\S1_total_effect_MR\harmonised_ebi-a-GCST006414.csv")["SNP"])
mat = piv.reset_index().rename(columns={"index": "rsid"})
mat = mat[mat["rsid"].isin(ref40)].copy()   # drop API alias rows (e.g. rs1157821370)
cat_cols = [c for c in mat.columns if c not in ("rsid", "total_hits")]
rows = []
for _, r in mat.iterrows():
    top = sorted(((c, int(r[c])) for c in cat_cols), key=lambda t: -t[1])
    top = [t for t in top if t[1] > 0]
    t1 = top[0] if top else ("—", 0)
    t2 = top[1] if len(top) > 1 else ("—", 0)
    rows.append({"rsid": r["rsid"], "total_hits": int(r["total_hits"]),
                 "top1": f"{t1[0]} ({t1[1]})", "top2": f"{t2[0]} ({t2[1]})"})
out = pd.DataFrame(rows).sort_values("total_hits", ascending=False)
out.to_csv(r"C:\Users\曹泽众\Documents\kimi\workspace\EJPC_修订产出\stats\phewas_snp_top2.csv", index=False)
print("exported phewas_category_summary.csv and phewas_snp_top2.csv (", len(out), "SNPs )")
