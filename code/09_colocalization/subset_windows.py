# -*- coding: utf-8 -*-
"""Subset FinnGen(GRCh38)/Nielsen(GRCh37) locus windows from ±1Mb to ±500kb
around lead SNPs, and emit the rsID list to fetch from OpenGWAS (ukb-b-5237)."""
import json
import pandas as pd

BASE = r"C:\Users\64915\Documents\kimi\workspace\integrated_v4\coloc_v42"
HALF = 500_000

# lead-SNP coordinates (verified previous session)
LEAD = {
    "FTO":           {"fg": 53767042, "nei": 53800954},   # rs1421085
    "TMEM18":        {"fg": 637498,   "nei": 637498},     # rs13387939
    "CYP1A1_CYP1A2": {"fg": 74735539, "nei": 75027880},   # rs2472297
    "AHR":           {"fg": 17244953, "nei": 17284577},   # rs4410790
    "PCMTD2":        {"fg": None,     "nei": 62891820},   # rs6062682 (absent in FG)
}

fg = pd.read_csv(f"{BASE}/fg_windows.tsv", sep="\t")
nei = pd.read_csv(f"{BASE}/nei_windows.tsv", sep="\t")

# PCMTD2 FG anchor: centre of the previously extracted ±1Mb window
pc = fg.loc[fg.locus == "PCMTD2", "pos"]
LEAD["PCMTD2"]["fg"] = int((pc.min() + pc.max()) // 2)
print("PCMTD2 FG anchor:", LEAD["PCMTD2"]["fg"], " (window", pc.min(), "-", pc.max(), ")")

fg_out, nei_out = [], []
for loc, c in LEAD.items():
    f = fg[(fg.locus == loc) & (fg.pos.between(c["fg"] - HALF, c["fg"] + HALF))]
    n = nei[(nei.locus == loc) & (nei.POS_GRCh37.between(c["nei"] - HALF, c["nei"] + HALF))]
    fg_out.append(f); nei_out.append(n)
    print(f"{loc:15s} chr_fg={sorted(f['#chrom'].unique())} chr_nei={sorted(n.CHR.unique())} "
          f"fg_rows={len(f):6d} nei_rows={len(n):6d}")

fg500 = pd.concat(fg_out); nei500 = pd.concat(nei_out)
fg500.to_csv(f"{BASE}/fg_windows_500kb.tsv", sep="\t", index=False)
nei500.to_csv(f"{BASE}/nei_windows_500kb.tsv", sep="\t", index=False)

rsids = sorted(set(fg500.rsids.dropna()))
with open(f"{BASE}/fetch_rsids.json", "w") as fh:
    json.dump(rsids, fh)
print("total fg rows:", len(fg500), "| total nei rows:", len(nei500), "| rsids to fetch:", len(rsids))
