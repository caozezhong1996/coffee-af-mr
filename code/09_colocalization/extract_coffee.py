# -*- coding: utf-8 -*-
"""Stream-parse ukb-b-5237.vcf.gz; extract rows for the 36,341 window rsIDs.
Output schema matches the old OpenGWAS API cache (id,trait,chr,position,rsid,ea,nea,eaf,beta,se,p,n)."""
import gzip, json
import pandas as pd

BASE = r"C:\Users\64915\Documents\kimi\workspace\integrated_v4\coloc_v42"
want = set(json.load(open(f"{BASE}/fetch_rsids.json")))
rows = []
bad = 0
with gzip.open(f"{BASE}/ukb-b-5237.vcf.gz", "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        c = line.split("\t", 9)
        rs = c[2]
        if rs not in want:
            continue
        try:
            es, se, lp, af, _id = c[9].rstrip("\n").split(":")
            rows.append((rs, c[0], int(c[1]), c[4], c[3], float(af),
                         float(es), float(se), 10.0 ** (-float(lp))))
        except Exception:
            bad += 1
print("matched rows:", len(rows), "| unparseable:", bad)

df = pd.DataFrame(rows, columns=["rsid", "chr", "position", "ea", "nea", "eaf", "beta", "se", "p"])
df.insert(0, "id", "ukb-b-5237")
df.insert(1, "trait", "Coffee intake")
df["n"] = 428860
df = df[["id", "trait", "chr", "position", "rsid", "ea", "nea", "eaf", "beta", "se", "p", "n"]]
df = df.drop_duplicates(subset=["id", "rsid"])
df.to_csv(f"{BASE}/coffee_vcf_associations.csv", index=False)
print("written:", len(df), "| coverage:", df.rsid.nunique(), "/", len(want))

fg = pd.read_csv(f"{BASE}/fg_windows_500kb.tsv", sep="\t")
got = set(df.rsid)
for loc in ["FTO", "TMEM18", "CYP1A1_CYP1A2", "AHR", "PCMTD2"]:
    w = set(fg.loc[fg.locus == loc, "rsids"].dropna())
    print(f"{loc:15s} fg={len(w):6d} coffee-matched={len(w & got):6d} ({100*len(w&got)/len(w):.0f}%)")
