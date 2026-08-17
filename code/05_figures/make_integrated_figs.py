# -*- coding: utf-8 -*-
import sys
from pathlib import Path
# Optional: configure CJK fonts via matplotlib rcParams if needed
# (font setup placeholder)
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

OUT = "integrated_paper"
import os; os.makedirs(OUT, exist_ok=True)

# ---------- Figure 1: integrated design ----------
fig, ax = plt.subplots(figsize=(10.5, 5.6))
ax.set_xlim(0, 10.5); ax.set_ylim(0, 6.1); ax.axis("off")
def box(x, y, w, h, text, fc, fs=8.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.12", fc=fc, ec="#333", lw=1.2))
    ax.text(x+w/2, y+h/2, text, ha="center", va="center", fontsize=fs, color="white", fontweight="bold", linespacing=1.3)
def arr(x1, y1, x2, y2, col="#333"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=14, color=col, lw=1.8))

box(3.9, 4.6, 2.7, 0.85, "Coffee / caffeine exposure", "#2c3e50", 9.5)
box(0.15, 2.6, 3.1, 1.35, "Component 1\nTotal-effect MR\nUKB coffee GWAS (40 SNP)\n→ Nielsen 2018 + FinnGen R11\n(116,473 AF cases)", "#1a5276")
box(3.7, 2.6, 3.1, 1.35, "Component 2\nNetwork MR (mechanism)\n6 prespecified mediators\nSBP · BMI · MVPA · CRP\nsleep · insomnia", "#7d3c98")
box(7.25, 2.6, 3.1, 1.35, "Component 3\nPharmacovigilance\nFAERS 11.9M serious reports\ncaffeine / citrate / taurine /\nenergy drinks (ROR, dedup, masking)", "#7d6608")
box(3.35, 0.25, 3.8, 1.15, "ATRIAL FIBRILLATION\ntotal effect · mechanism · real-world safety", "#7b241c", 9)
arr(5.25, 4.6, 1.7, 3.95); arr(5.25, 4.6, 5.25, 3.95); arr(5.25, 4.6, 8.8, 3.95)
arr(1.7, 2.6, 4.2, 1.4); arr(5.25, 2.6, 5.25, 1.4); arr(8.8, 2.6, 6.4, 1.4)
ax.text(5.25, 5.9, "Integrated design: triangulation of total effect, mediating pathways, and real-world safety",
        ha="center", fontsize=11.5, fontweight="bold")
fig.savefig(f"{OUT}/Fig1_overview_300dpi.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig1_overview_300dpi.tiff", dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------- Figure 3: total-effect forest across specifications ----------
rows = [
 ("Nielsen 2018 — univariable IVW", 1.22, 0.97, 1.53),
 ("FinnGen R11 — univariable IVW", 1.63, 1.13, 2.36),
 ("Pooled FE (Nielsen + R11)", 1.32, 1.09, 1.61),
 ("Pooled RE (Nielsen + R11)", 1.36, 1.03, 1.80),
 ("Pooled MR-PRESSO-corrected FE", 1.19, 1.02, 1.40),
 ("Pooled MR-PRESSO-corrected RE", 1.22, 0.94, 1.58),
 ("Nielsen — MVMR (BMI, smoking, alcohol)", 1.01, 0.81, 1.26),
 ("FinnGen R11 — MVMR (BMI, smoking, alcohol)", 1.16, 0.85, 1.59),
 ("Pooled MVMR", 1.06, 0.88, 1.27),
]
fig, ax = plt.subplots(figsize=(8.8, 4.9))
ys = range(len(rows))[::-1]
for y, (lab, orr, lo, hi) in zip(ys, rows):
    sig = lo > 1 or hi < 1
    col = "#c0392b" if sig else "#1a5276" if hi < 1 else "#555555"
    ax.plot([lo, hi], [y, y], color=col, lw=2)
    ax.scatter([orr], [y], s=60, color=col, zorder=3)
    ax.text(2.62, y, f"{orr:.2f} ({lo:.2f}–{hi:.2f})" + ("*" if sig else ""), va="center", fontsize=8.5)
ax.axvline(1, color="k", lw=0.8, ls="--")
ax.set_yticks(list(ys)); ax.set_yticklabels([r[0] for r in rows], fontsize=9.5)
ax.set_xscale("log"); ax.set_xlim(0.7, 2.6)
ax.set_xticks([0.8, 1.0, 1.25, 1.5, 2.0, 2.5]); ax.set_xticklabels(["0.8","1.0","1.25","1.5","2.0","2.5"])
ax.set_xlabel("OR for atrial fibrillation per +1 cup/day of genetically predicted coffee intake (log scale)")
ax.set_title("Total-effect MR estimates across datasets and pleiotropy-control specifications", fontsize=11)
fig.savefig(f"{OUT}/Fig3_totaleffect_forest_300dpi.png", dpi=300, bbox_inches="tight")
fig.savefig(f"{OUT}/Fig3_totaleffect_forest_300dpi.tiff", dpi=300, bbox_inches="tight")
plt.close(fig)
print("done")
