# -*- coding: utf-8 -*-
import json, sys
from pathlib import Path
# Optional: configure CJK fonts via matplotlib rcParams if needed
# (font setup placeholder)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = "mediation_mr"
S = json.load(open(f"{OUT}/summary_stats.json"))
MEDS = ["sbp","bmi","mvpa","crp","sleep","insomnia"]
LBL = dict(sbp="Systolic blood pressure", bmi="Body mass index", mvpa="Moderate-to-vigorous\nphysical activity",
           crp="C-reactive protein", sleep="Sleep duration", insomnia="Insomnia")
UNIT = dict(sbp="mmHg", bmi="kg/m²", mvpa="SD", crp="log mg/L", sleep="h", insomnia="log OR")

# ---------------- Figure 2: step1 forest ----------------
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ys = np.arange(len(MEDS))[::-1]
for y, k in zip(ys, MEDS):
    d = S["step1"][k]
    b, se, p = d["ivw_b"], d["ivw_se"], d["ivw_p"]
    lo, hi = b-1.96*se, b+1.96*se
    col = "#c0392b" if (lo > 0) else ("#1a5276" if hi < 0 else "#555555")
    ax.plot([lo, hi], [y, y], color=col, lw=2)
    ax.scatter([b], [y], s=60, color=col, zorder=3)
    sig = "*" if p < 0.05 else ""
    ax.text(1.25, y, f"β={b:+.3f} ({lo:+.3f} to {hi:+.3f}){sig}", va="center", fontsize=8.5)
ax.axvline(0, color="k", lw=0.8, ls="--")
ax.set_yticks(ys); ax.set_yticklabels([f"{LBL[k]}  [{UNIT[k]}]" for k in MEDS], fontsize=9.5)
ax.set_xlim(-2.6, 1.3)
ax.set_xlabel("Effect of genetically predicted coffee intake (per +1 cup/day) on mediator, β (95% CI)")
ax.set_title("Step 1 MR: coffee consumption → candidate mediators (IVW)", fontsize=11)
for ext in ["png","tiff"]:
    fig.savefig(f"{OUT}/Fig2_step1_forest_300dpi.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------------- Figure 3: step2 forest ----------------
fig, axes = plt.subplots(1, 2, figsize=(10, 4.8), sharey=True)
for ax, oc, ttl in zip(axes, ["nielsen","r11"], ["AF (Nielsen 2018, 60,620 cases)", "AF/flutter (FinnGen R11, 45,766 cases)"]):
    for y, k in zip(ys, MEDS):
        d = S["step2"][k][oc]
        b, se, p = d["ivw_b"], d["ivw_se"], d["ivw_p"]
        lo, hi = b-1.96*se, b+1.96*se
        col = "#c0392b" if (lo > 0) else ("#1a5276" if hi < 0 else "#555555")
        ax.plot([lo, hi], [y, y], color=col, lw=2)
        ax.scatter([b], [y], s=60, color=col, zorder=3)
        sig = "*" if p < 0.05 else ""
        ax.text(0.85, y, f"{b:+.3f} ({lo:+.3f},{hi:+.3f}){sig}", va="center", fontsize=8)
    ax.axvline(0, color="k", lw=0.8, ls="--")
    ax.set_xlim(-0.55, 1.0)
    ax.set_xlabel("Causal effect on AF, log OR (95% CI)")
    ax.set_title(ttl, fontsize=10)
axes[0].set_yticks(ys); axes[0].set_yticklabels([LBL[k] for k in MEDS], fontsize=9.5)
fig.suptitle("Step 2 MR: mediators → atrial fibrillation (IVW)", fontsize=11)
fig.tight_layout(rect=[0,0,1,0.94])
for ext in ["png","tiff"]:
    fig.savefig(f"{OUT}/Fig3_step2_forest_300dpi.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------------- Figure 4: indirect effects ----------------
ind = [x for x in S["indirect"] if x["step1_method"]=="ivw"]
fig, ax = plt.subplots(figsize=(7.6, 5.2))
labels = []
ys2 = []
y = 0
pos = []
for k in MEDS:
    for oc, ol in [("nielsen","Nielsen"), ("r11","FinnGen")]:
        x = [i for i in ind if i["mediator"]==k and i["outcome"]==oc][0]
        b = float(x["indirect"]); se = float(x["indirect_se"]); p = float(x["p"])
        lo, hi = b-1.96*se, b+1.96*se
        pos.append((y, b, lo, hi, p, oc, k))
        labels.append((y, k, oc))
        y -= 1
    y -= 0.6
for yy, b, lo, hi, p, oc, k in pos:
    col = "#c0392b" if (lo > 0) else ("#1a5276" if hi < 0 else "#555555")
    marker = "o" if oc == "nielsen" else "s"
    ax.plot([lo, hi], [yy, yy], color=col, lw=2)
    ax.scatter([b], [yy], s=55, color=col, marker=marker, zorder=3)
    sig = "*" if p < 0.05 else ""
    ax.text(0.62, yy, f"{b:+.3f} ({lo:+.3f},{hi:+.3f}){sig}", va="center", fontsize=8)
ax.axvline(0, color="k", lw=0.8, ls="--")
ax.set_yticks([l[0] for l in labels])
ax.set_yticklabels([f"{LBL[k].replace(chr(10),' ')}  —  {'Nielsen' if oc=='nielsen' else 'FinnGen R11'}" for _, k, oc in labels], fontsize=8.5)
ax.set_xlim(-0.55, 0.68)
ax.set_xlabel("Indirect (mediated) effect of coffee on AF, log OR (95% CI, Sobel)")
ax.set_title("Indirect effects of coffee on atrial fibrillation via each mediator\n(product-of-coefficients, IVW × IVW)", fontsize=11)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0],[0], marker="o", color="w", markerfacecolor="k", label="Nielsen 2018", markersize=7),
                   Line2D([0],[0], marker="s", color="w", markerfacecolor="k", label="FinnGen R11", markersize=7)],
          loc="lower right", fontsize=8.5)
for ext in ["png","tiff"]:
    fig.savefig(f"{OUT}/Fig4_indirect_forest_300dpi.{ext}", dpi=300, bbox_inches="tight")
plt.close(fig)
print("Fig2-4 done")
