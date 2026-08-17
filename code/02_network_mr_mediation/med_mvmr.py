# -*- coding: utf-8 -*-
"""Manual multivariable IVW: coffee + mediator -> AF, using local harmonised data."""
import gzip, math, os
import numpy as np
import pandas as pd

OUT = "mediation_mr"
MEDS = ["sleep", "insomnia", "sbp", "mvpa", "crp", "bmi"]
PAL = {("A","T"), ("T","A"), ("G","C"), ("C","G")}

# ---------- load sources ----------
inst = pd.read_csv("mr_results/instruments.csv").to_dict("records")
coffee_ref = {r["SNP"]: (r["effect_allele.exposure"].upper(), r["other_allele.exposure"].upper(),
                      r["beta.exposure"], r["se.exposure"], r.get("eaf.exposure", np.nan))
              for r in inst}

# coffee betas at union snps from API
union = set(l.strip() for l in open(os.path.join(OUT, "union_snps.txt")))
cr = pd.read_csv(os.path.join(OUT, "coffee_union_raw.csv")).to_dict("records")
cgwas = {r["rsid"]: (str(r["ea"]).upper(), str(r["nea"]).upper(), r["beta"], r["se"], r["eaf"]) for r in cr}
print("coffee gwas coverage:", len(cgwas), "/", len(union))

niel = pd.read_csv(os.path.join(OUT, "outcome_nielsen_union_raw.csv")).to_dict("records")
niel_ref = {r["rsid"]: (str(r["ea"]).upper(), str(r["nea"]).upper(), r["beta"], r["se"], r["eaf"]) for r in niel}

r11 = pd.read_csv(os.path.join(OUT, "r11_union_hits.tsv"), sep="\t")
r11_ref = {}
for r in r11.itertuples():
    hits = [s for s in str(r.rsids).split(",") if s in union]
    if hits:
        r11_ref[hits[0]] = (str(r.alt).upper(), str(r.ref).upper(), r.beta, r.sebeta, r.af_alt)

def align(rec, ea_ref, oa_ref):
    """rec = (EA, OA, beta, se, eaf). Align beta to ea_ref; return (beta, se) or None."""
    ea, oa, b, se, eaf = rec
    if ea == ea_ref and oa == oa_ref:
        return b, se, eaf
    if ea == oa_ref and oa == ea_ref:
        return -b, se, (1 - eaf if not (eaf is None or (isinstance(eaf, float) and math.isnan(eaf))) else np.nan)
    return None

def is_pal_amb(ea, oa, eaf):
    return (ea, oa) in PAL and (eaf is not None and not (isinstance(eaf, float) and math.isnan(eaf)) and 0.42 < eaf < 0.58)

def mvmr_ivw(df):
    X = df[["bx_coffee", "bx_med"]].values
    y = df["by"].values
    w = 1.0 / df["seY"].values ** 2
    XtW = X.T * w
    V = np.linalg.inv(XtW @ X)
    beta = V @ XtW @ y
    resid = y - X @ beta
    dfree = len(y) - 2
    sigma = max(1.0, math.sqrt((w * resid ** 2).sum() / dfree))
    se = np.sqrt(np.diag(V)) * sigma
    from scipy import stats
    z = beta / se
    p = 2 * stats.norm.sf(np.abs(z))
    return beta, se, p, sigma

rows = []
detail_dir = os.path.join(OUT, "mvmr_data")
os.makedirs(detail_dir, exist_ok=True)
for k in MEDS:
    keep = pd.read_csv(os.path.join(OUT, f"mvmr_snps_{k}.csv")).to_dict("records")
    iv = pd.read_csv(os.path.join(OUT, f"iv_{k}.csv")).to_dict("records")
    iv_ref = {r["SNP"]: (r["effect_allele.exposure"].upper(), r["other_allele.exposure"].upper(),
                      r["beta.exposure"], r["se.exposure"], r.get("eaf.exposure", np.nan))
              for r in iv}
    s1 = pd.read_csv(os.path.join(OUT, f"step1_harmonised_{k}.csv")).to_dict("records")
    s1_med = {r["SNP"]: (r["beta.outcome"], r["se.outcome"]) for r in s1 if r["mr_keep"]}

    recs = []
    for r in keep:
        snp = r["SNP"]
        if r["src"] == "coffee":
            if snp not in coffee_ref or snp not in s1_med: continue
            ea, oa, bc, sec, eafc = coffee_ref[snp]
            if is_pal_amb(ea, oa, eafc): continue
            bm, sem = s1_med[snp]
            row = dict(SNP=snp, EA=ea, OA=oa, bx_coffee=bc, se_coffee=sec, bx_med=bm, se_med=sem, eaf=eafc)
        else:
            if snp not in iv_ref or snp not in cgwas: continue
            ea, oa, bm, sem, eafm = iv_ref[snp]
            if is_pal_amb(ea, oa, eafm): continue
            a = align(cgwas[snp], ea, oa)
            if a is None: continue
            row = dict(SNP=snp, EA=ea, OA=oa, bx_coffee=a[0], se_coffee=a[1], bx_med=bm, se_med=sem, eaf=eafm)
        recs.append(row)
    base = pd.DataFrame(recs)
    for oc, oref in [("nielsen", niel_ref), ("r11", r11_ref)]:
        out_rows = []
        for r in base.to_dict("records"):
            if r["SNP"] not in oref: continue
            a = align(oref[r["SNP"]], r["EA"], r["OA"])
            if a is None: continue
            out_rows.append(dict(SNP=r["SNP"], bx_coffee=r["bx_coffee"], bx_med=r["bx_med"], by=a[0], seY=a[1]))
        df = pd.DataFrame(out_rows).dropna()
        if len(df) < 10:
            print(k, oc, "too few:", len(df)); continue
        df.to_csv(os.path.join(detail_dir, f"mvmr_{k}_{oc}.csv"), index=False)
        beta, se, p, sigma = mvmr_ivw(df)
        rows.append(dict(mediator=k, outcome=oc, nsnp=len(df), sigma=round(sigma, 3),
                         coffee_b=beta[0], coffee_se=se[0], coffee_p=p[0],
                         med_b=beta[1], med_se=se[1], med_p=p[1]))

res = pd.DataFrame(rows)
for c in ["coffee_or", "med_or"]: pass
res["coffee_or"] = np.exp(res.coffee_b).round(3)
res["coffee_or_lo"] = np.exp(res.coffee_b - 1.96*res.coffee_se).round(3)
res["coffee_or_hi"] = np.exp(res.coffee_b + 1.96*res.coffee_se).round(3)
res["med_or"] = np.exp(res.med_b).round(3)
res.to_csv(os.path.join(OUT, "mvmr_direct_effects.csv"), index=False)
pd.set_option("display.width", 200)
print(res[["mediator","outcome","nsnp","sigma","coffee_b","coffee_se","coffee_p","coffee_or","med_b","med_p"]].to_string(index=False))
print("SAVED mvmr_direct_effects.csv")
