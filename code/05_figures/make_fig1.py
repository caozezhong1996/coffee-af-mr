# -*- coding: utf-8 -*-
import sys, json
from pathlib import Path
# Optional: configure CJK fonts via matplotlib rcParams if needed
# (font setup placeholder)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

S = json.load(open("mediation_mr/summary_stats.json"))
fig, ax = plt.subplots(figsize=(9.6, 6.0))
ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")

def box(x, y, w, h, text, fc, fs=9.5, tc="white", bold=True):
    b = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.12", fc=fc, ec="#333333", lw=1.2)
    ax.add_patch(b)
    ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fs, color=tc,
            fontweight="bold" if bold else "normal", linespacing=1.35)

def arrow(x1, y1, x2, y2, label="", col="#333333", ls="-", lx=None, ly=None, fs=8.5):
    a = FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=14, color=col, lw=1.8, linestyle=ls)
    ax.add_patch(a)
    if label:
        ax.text(lx if lx is not None else (x1+x2)/2, ly if ly is not None else (y1+y2)/2, label,
                fontsize=fs, color=col, ha="center", va="center",
                bbox=dict(fc="white", ec="none", alpha=0.9, pad=1.2))

# main nodes
box(0.2, 4.6, 2.2, 1.1, "Habitual coffee\nconsumption\n(UKB, n=428,860;\n40 SNPs)", "#2c3e50", fs=8.5)
box(7.6, 4.6, 2.2, 1.1, "Atrial\nfibrillation\nNielsen 2018 +\nFinnGen R11", "#7b241c", fs=8.5)

# mediator nodes
med_y = dict(sbp=8.4, bmi=6.7, mvpa=5.0, crp=3.3, sleep=1.6, insomnia=-0.1)
lbl = dict(sbp="Systolic BP", bmi="BMI", mvpa="MVPA", crp="CRP", sleep="Sleep duration", insomnia="Insomnia")
col = dict(sbp="#1a5276", bmi="#922b21", mvpa="#117864", crp="#7d6608", sleep="#5d6d7e", insomnia="#5d6d7e")
for k, y in med_y.items():
    yy = max(y, 0.35)
    box(4.0, yy, 2.0, 0.85, lbl[k], col[k], fs=9)

b1 = dict(sbp=-1.934, bmi=+0.804, mvpa=-0.146, crp=-0.123, sleep=-0.013, insomnia=-0.418)
b2n = {k: S["step2"][k]["nielsen"]["ivw_b"] for k in med_y}
p1 = {k: S["step1"][k]["ivw_p"] for k in med_y}
p2 = {k: S["step2"][k]["nielsen"]["ivw_p"] for k in med_y}

for k, y in med_y.items():
    yy = max(y, 0.35) + 0.42
    c1 = "#1a5276" if b1[k] < 0 and p1[k] < 0.05 else ("#c0392b" if p1[k] < 0.05 else "#808b96")
    arrow(2.4, 5.15, 4.0, yy, col=c1, ls="-" if p1[k] < 0.05 else "--")
    ax.text(3.0, (5.15+yy)/2, f"β1={b1[k]:+.2f}" + ("*" if p1[k] < 0.05 else ""), fontsize=8, color=c1,
            ha="center", va="center", bbox=dict(fc="white", ec="none", alpha=0.9, pad=1.0))
    c2 = "#c0392b" if b2n[k] > 0 and p2[k] < 0.05 else ("#1a5276" if p2[k] < 0.05 else "#808b96")
    arrow(6.0, yy, 7.6, 5.15, col=c2, ls="-" if p2[k] < 0.05 else "--")
    ax.text(7.0, (5.15+yy)/2, f"β2={b2n[k]:+.3f}" + ("*" if p2[k] < 0.05 else ""), fontsize=8, color=c2,
            ha="center", va="center", bbox=dict(fc="white", ec="none", alpha=0.9, pad=1.0))

# total effect
arrow(2.4, 4.75, 7.6, 4.75, col="#555555", ls=":")
ax.text(0.3, 9.3, "Total effect: OR 1.22 (0.97–1.53) Nielsen;\nOR 1.63 (1.13–2.35) FinnGen R11;\nOR≈1.01 after confounder MVMR",
        fontsize=8, color="#555555", ha="left", va="top",
        bbox=dict(fc="#f8f9f9", ec="#cccccc", pad=4))

ax.text(5.0, 9.7, "Two-step network Mendelian randomization: design and key path coefficients",
        fontsize=11.5, ha="center", fontweight="bold")
ax.text(5.0, 0.05, "Solid arrows: P<0.05; dashed: P≥0.05.  β1: per +1 cup/day; β2: log OR per mediator unit (Nielsen).  *P<0.05",
        fontsize=7.5, ha="center", color="#666666")
for ext in ["png","tiff"]:
    fig.savefig("mediation_mr/Fig1_design_300dpi."+ext, dpi=300, bbox_inches="tight")
plt.close(fig)
print("Fig1 done")
