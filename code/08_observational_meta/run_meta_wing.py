# -*- coding: utf-8 -*-
"""Updated dose-response meta-analysis: coffee intake (cups/day) and incident AF.
Primary: per-cup linear trend — Greenland-Longnecker GLS for cohorts with per-category
cases (COSM, SMC, DCH, PHS) + directly reported per-cup HR (Kim 2021 UK Biobank),
pooled with DerSimonian-Laird random effects; leave-one-out sensitivity.
Secondary: highest-vs-lowest category pooling (broader study set incl. caffeine-based).
All numbers extracted from PMC full texts / abstracts (see extractions_meta_wing.csv)."""
import io, sys, math, json
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
import numpy as np
import pandas as pd

Z = 1.959964

def se_from_ci(lo, hi):
    return (math.log(hi) - math.log(lo)) / (2 * Z)

def gl_trend(doses, rrs, los, his, cases_ref):
    """Greenland-Longnecker GLS linear trend per 1 unit dose (through origin)."""
    d = np.array(doses, float)
    y = np.array([0.0] + [math.log(r) for r in rrs])
    s0 = math.sqrt(1.0 / cases_ref)  # var of log rate in reference category
    s = [s0] + [se_from_ci(lo, hi) for lo, hi in zip(los, his)]
    S = np.diag(np.square(s))
    S += s0 ** 2  # shared reference covariance (off-diagonal)
    S[0, 0] = s0 ** 2  # ref cell diagonal
    Si = np.linalg.inv(S)
    dX = d.reshape(-1, 1)
    beta = ((dX.T @ Si @ y) / (dX.T @ Si @ dX)).item()
    var = (1.0 / (dX.T @ Si @ dX)).item()
    return beta, math.sqrt(var)

studies = {}

# --- Larsson 2015 COSM (men); doses: <2→1, 2-<3→2.5, 3-<4→3.5, 4-<5→4.5, ≥5→6
studies["Larsson 2015 (COSM, men)"] = gl_trend(
    [1, 2.5, 3.5, 4.5, 6],
    [0.99, 0.99, 0.97, 1.08],
    [0.90, 0.89, 0.87, 0.98],
    [1.09, 1.09, 1.07, 1.20], 742)
# --- Larsson 2015 SMC (women)
studies["Larsson 2015 (SMC, women)"] = gl_trend(
    [1, 2.5, 3.5, 4.5, 6],
    [0.95, 0.98, 0.98, 0.88],
    [0.85, 0.87, 0.86, 0.76],
    [1.07, 1.10, 1.12, 1.02], 505)
# --- Mostofsky 2016 Danish DCH; doses 0,.5,1,2.5,4.5,6.5,8
studies["Mostofsky 2016 (DCH)"] = gl_trend(
    [0, 0.5, 1, 2.5, 4.5, 6.5, 8],
    [0.93, 0.88, 0.86, 0.84, 0.79, 0.79],
    [0.74, 0.71, 0.71, 0.69, 0.64, 0.63],
    [1.15, 1.10, 1.04, 1.02, 0.98, 1.00], 124)
# --- Physicians' Health Study 2019 (Model 3)
studies["Bodar 2019 (PHS)"] = gl_trend(
    [0, 0.14, 0.43, 0.79, 1, 2.5, 4.5],
    [0.84, 1.05, 0.88, 0.87, 0.85, 0.93],
    [0.69, 0.86, 0.69, 0.75, 0.75, 0.77],
    [1.02, 1.29, 1.13, 1.00, 0.97, 1.12], 472)
# --- Kim 2021 UK Biobank: directly reported per-cup HR
studies["Kim 2021 (UK Biobank)"] = (math.log(0.97), se_from_ci(0.96, 0.98))

rows = []
for k, (b, s) in studies.items():
    rows.append({"study": k, "logHR_per_cup": b, "se": s,
                 "HR": math.exp(b), "lo": math.exp(b - Z * s), "hi": math.exp(b + Z * s)})
df = pd.DataFrame(rows)

# DerSimonian-Laird random effects
w = 1 / df["se"] ** 2
b_fe = float((w * df["logHR_per_cup"]).sum() / w.sum())
Q = float((w * (df["logHR_per_cup"] - b_fe) ** 2).sum())
k = len(df)
C = w.sum() - (w ** 2).sum() / w.sum()
tau2 = max(0.0, (Q - (k - 1)) / C)
wr = 1 / (df["se"] ** 2 + tau2)
b_re = float((wr * df["logHR_per_cup"]).sum() / wr.sum())
se_re = math.sqrt(1 / wr.sum())
I2 = max(0.0, (Q - (k - 1)) / Q) * 100
res = {"k": k, "Q": Q, "Q_df": k - 1, "I2_pct": I2, "tau2": tau2,
       "HR_per_cup_RE": math.exp(b_re),
       "lo": math.exp(b_re - Z * se_re), "hi": math.exp(b_re + Z * se_re),
       "P": 2 * (1 - 0.5 * (1 + math.erf(abs(b_re / se_re) / math.sqrt(2))))}

# leave-one-out
loo = []
for i in range(k):
    sub = df.drop(i)
    w2 = 1 / sub["se"] ** 2
    bf = float((w2 * sub["logHR_per_cup"]).sum() / w2.sum())
    Q2 = float((w2 * (sub["logHR_per_cup"] - bf) ** 2).sum())
    C2 = w2.sum() - (w2 ** 2).sum() / w2.sum()
    t2 = max(0.0, (Q2 - (k - 2)) / C2)
    wr2 = 1 / (sub["se"] ** 2 + t2)
    br = float((wr2 * sub["logHR_per_cup"]).sum() / wr2.sum())
    sr = math.sqrt(1 / wr2.sum())
    loo.append({"omitted": df["study"][i], "HR": math.exp(br),
                "lo": math.exp(br - Z * sr), "hi": math.exp(br + Z * sr)})
loo = pd.DataFrame(loo)

# ---- Secondary: highest vs lowest category ----
hl = pd.DataFrame([
    # study, RR, lo, hi, exposure note
    ["Wilhelmsen 2001 (Gothenburg)", 1.09, 0.87, 1.38, "≥5 vs 0 cups/d"],
    ["Conen 2010 (WHS)", 1.03, 0.79, 1.35, "≥4 vs 0 cups/d"],
    ["Klatsky 2011 (Kaiser)", 0.81, 0.69, 0.96, "≥4 vs 0 cups/d"],
    ["Larsson 2015 (COSM, men)", 1.08, 0.98, 1.20, "≥5 vs <2 cups/d"],
    ["Larsson 2015 (SMC, women)", 0.88, 0.76, 1.02, "≥5 vs <2 cups/d"],
    ["Mostofsky 2016 (DCH)", 0.79, 0.63, 1.00, ">7 vs 0 cups/d"],
    ["Bodar 2019 (PHS)", 0.93, 0.77, 1.12, "4+ vs almost never"],
    ["Bazal 2021 (SUN+PREDIMED)", 0.79, 0.49, 1.28, ">1 vs ≤3/mo cups"],
    ["Sehrawat 2023 (MESA)", 1.40, 1.07, 1.84, "≥1 cup/wk vs none (SE from P=0.015)"],
    ["Mattina? 2018 (Italian cohort)", 0.249, 0.161, 0.458, "3rd vs 1st caffeine tertile"],
], columns=["study", "RR", "lo", "hi", "contrast"])
hl["logRR"] = np.log(hl["RR"]); hl["se"] = (np.log(hl["hi"]) - np.log(hl["lo"])) / (2 * Z)
wh = 1 / hl["se"] ** 2
bf = float((wh * hl["logRR"]).sum() / wh.sum())
Qh = float((wh * (hl["logRR"] - bf) ** 2).sum())
Ch = wh.sum() - (wh ** 2).sum() / wh.sum()
t2h = max(0.0, (Qh - (len(hl) - 1)) / Ch)
wrh = 1 / (hl["se"] ** 2 + t2h)
brh = float((wrh * hl["logRR"]).sum() / wrh.sum()); srh = math.sqrt(1 / wrh.sum())
hlpool = {"k": len(hl), "Q": Qh, "I2_pct": max(0.0, (Qh - (len(hl) - 1)) / Qh) * 100,
          "RR_RE": math.exp(brh), "lo": math.exp(brh - Z * srh), "hi": math.exp(brh + Z * srh)}
# without the two caffeine-based / small outliers
hl2 = hl.drop([8, 9])
wh2 = 1 / hl2["se"] ** 2
bf2 = float((wh2 * hl2["logRR"]).sum() / wh2.sum())
Qh2 = float((wh2 * (hl2["logRR"] - bf2) ** 2).sum())
Ch2 = wh2.sum() - (wh2 ** 2).sum() / wh2.sum()
t2h2 = max(0.0, (Qh2 - (len(hl2) - 1)) / Ch2)
wrh2 = 1 / (hl2["se"] ** 2 + t2h2)
brh2 = float((wrh2 * hl2["logRR"]).sum() / wrh2.sum()); srh2 = math.sqrt(1 / wrh2.sum())
hlpool2 = {"k": len(hl2), "Q": Qh2, "I2_pct": max(0.0, (Qh2 - (len(hl2) - 1)) / max(Qh2, 1e-9)) * 100,
           "RR_RE": math.exp(brh2), "lo": math.exp(brh2 - Z * srh2), "hi": math.exp(brh2 + Z * srh2)}

df.to_csv("percup_trends.csv", index=False)
hl.to_csv("high_vs_low.csv", index=False)
loo.to_csv("loo_percup.csv", index=False)
json.dump({"percup_RE": res, "highlow_RE": hlpool, "highlow_RE_no_outliers": hlpool2},
          open("meta_wing_results.json", "w"), indent=1)

print("=== PRIMARY: per-cup linear trend (random effects) ===")
print(df.assign(HR=df.HR.round(3), lo=df.lo.round(3), hi=df.hi.round(3)).to_string(index=False))
print(f"Pooled HR per cup/day = {res['HR_per_cup_RE']:.3f} ({res['lo']:.3f}-{res['hi']:.3f}), "
      f"P={res['P']:.4f}, I2={I2:.0f}%, Q P df={k-1}")
print("\n=== Leave-one-out ===")
print(loo.round(3).to_string(index=False))
print(f"\n=== High-vs-low (all {hlpool['k']}): RR {hlpool['RR_RE']:.3f} ({hlpool['lo']:.3f}-{hlpool['hi']:.3f}), I2={hlpool['I2_pct']:.0f}%")
print(f"=== High-vs-low (excl. MESA & Italian): RR {hlpool2['RR_RE']:.3f} ({hlpool2['lo']:.3f}-{hlpool2['hi']:.3f}), I2={hlpool2['I2_pct']:.0f}%")
