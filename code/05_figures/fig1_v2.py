from pathlib import Path
import sys
# Optional: configure CJK fonts via matplotlib rcParams if needed
# (font setup placeholder)

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

NAVY = "#1f3a5f"; BLUE = "#2f6db3"; GREEN = "#1e7d46"; RED = "#b03a2e"; GRAY = "#6b7280"; LIGHT = "#eef3f8"

def box(ax, x, y, w, h, text, fc="white", ec=NAVY, fs=8.2, lw=1.4, tc="black"):
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012,rounding_size=0.015",
                       linewidth=lw, edgecolor=ec, facecolor=fc, mutation_aspect=1)
    ax.add_patch(p)
    ax.text(x + w/2, y + h/2, text, ha="center", va="center", fontsize=fs, color=tc, linespacing=1.45)

def arrow(ax, xy1, xy2, color=NAVY, lw=1.6, style="-|>", ls="-", rad=0.0):
    a = FancyArrowPatch(xy1, xy2, arrowstyle=style, mutation_scale=14,
                        linewidth=lw, color=color, linestyle=ls,
                        connectionstyle=f"arc3,rad={rad}", shrinkA=2, shrinkB=2)
    ax.add_patch(a)

fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.6, 7.2), gridspec_kw={"width_ratios": [1.05, 1.25]})
for ax in (axA, axB):
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

# ---------------- Panel A: triangulation framework ----------------
axA.text(0.02, 0.99, "A", fontsize=15, fontweight="bold", ha="left", va="top")
box(axA, 0.13, 0.885, 0.74, 0.085,
    "Does habitual coffee / caffeine intake\ncause atrial fibrillation?",
    fc=LIGHT, fs=9.6, lw=1.8)

pillars = [
    (0.005, "Total-effect MR\n(genetic causal inference)",
     "Pooled IVW OR 1.32 (1.09-1.61)\n-> MR-PRESSO 1.19 -> MVMR 1.06\nTOST-equivalent (P = 0.037)\nColoc: FTO/TMEM18 shared;\nCYP1A1 distinct (linkage)",
     "Key bias: horizontal pleiotropy"),
    (0.345, "Network mediation MR\n(mechanistic decomposition)",
     "Protective via SBP:\n-0.032 / -0.109 log OR\nDeleterious via BMI*:\n+0.291 / +0.477 log OR\nAdj. SBP unmasks\nOR 1.44-1.90",
     "Key bias: instrument proxy\n(Steiger: BMI direction n.s.)"),
    (0.685, "FAERS pharmacovigilance\n(real-world safety)",
     "11.9M serious reports\nCaffeine ROR 0.81 (0.68-0.96)\nNo signal after deduplication\nand masking correction",
     "Key bias: reporting bias /\nunder-reporting"),
]
for x0, title, body, bias in pillars:
    box(axA, x0, 0.635, 0.31, 0.085, title, fc="white", fs=8.4, lw=1.5)
    box(axA, x0, 0.30, 0.31, 0.285, body, fc="#fbfcfe", ec=BLUE, fs=7.6, lw=1.1)
    axA.text(x0 + 0.155, 0.265, bias, ha="center", va="top", fontsize=6.9, color=GRAY, style="italic", linespacing=1.3)
    arrow(axA, (x0 + 0.155, 0.63), (x0 + 0.155, 0.59), color=BLUE, lw=1.2)

for x0 in (0.005, 0.345, 0.685):
    arrow(axA, (0.5, 0.88), (x0 + 0.155, 0.725), color=NAVY, lw=1.3)
    arrow(axA, (x0 + 0.155, 0.225), (0.5, 0.15), color=NAVY, lw=1.3)

box(axA, 0.09, 0.02, 0.82, 0.125,
    "Convergent interpretation\nCustomary coffee intake is unlikely to be a clinically important cause of AF.\nThe fragile net association reflects offsetting pathways, not biological indifference.\nRoutine caffeine restriction for AF prevention is not supported.",
    fc="#eaf4ec", ec=GREEN, fs=8.4, lw=1.6)

# ---------------- Panel B: network MR path diagram ----------------
axB.text(0.0, 0.99, "B", fontsize=15, fontweight="bold", ha="left", va="top")

# main nodes
box(axB, 0.02, 0.50, 0.20, 0.11, "Genetically predicted\ncoffee intake\n(cups/day; 40 SNPs)", fc=LIGHT, fs=8.6, lw=1.6)
box(axB, 0.78, 0.50, 0.205, 0.11, "Atrial fibrillation\nNielsen 2018 (60,620)\nFinnGen R11 (55,853)", fc=LIGHT, fs=8.6, lw=1.6)

# total effect arrow
arrow(axB, (0.225, 0.555), (0.775, 0.555), color=GRAY, lw=1.8)
axB.text(0.5, 0.575, "Total effect (fragile): pooled IVW OR 1.32 -> MVMR 1.06\n(null-equivalent within OR 0.80-1.25; TOST P = 0.037)",
         ha="center", va="bottom", fontsize=7.6, color=GRAY)

# SBP pathway (top)
box(axB, 0.415, 0.72, 0.17, 0.085, "Systolic blood\npressure (mmHg)", fc="#eaf4ec", ec=GREEN, fs=8.4, lw=1.5)
arrow(axB, (0.13, 0.615), (0.455, 0.715), color=GREEN, lw=1.7)
arrow(axB, (0.545, 0.715), (0.87, 0.615), color=GREEN, lw=1.7)
axB.text(0.24, 0.75, "beta1 = -1.93\n(P = 0.027)", ha="center", va="center", fontsize=7.6, color=GREEN)
axB.text(0.72, 0.73, "beta2 = +0.016 / +0.056\n(P = 2e-14 / 3e-85)", ha="center", va="center", fontsize=7.6, color=GREEN)
axB.text(0.46, 0.825, "Indirect effect -0.032 / -0.109 (protective)\nP = 0.033 / 0.028", ha="right", va="bottom", fontsize=8.2, color=GREEN, fontweight="bold")

# null mediators (top right)
box(axB, 0.63, 0.90, 0.20, 0.075, "Sleep / insomnia /\nCRP / MVPA", fc="#f3f4f6", ec=GRAY, fs=7.8, lw=1.2)
arrow(axB, (0.17, 0.615), (0.65, 0.935), color=GRAY, lw=1.0, ls="--")
arrow(axB, (0.83, 0.935), (0.90, 0.615), color=GRAY, lw=1.0, ls="--")
axB.text(0.945, 0.925, "No mediation\n(P > 0.48)", ha="center", va="center", fontsize=7.2, color=GRAY)

# BMI pathway (bottom)
box(axB, 0.415, 0.24, 0.17, 0.085, "Body mass\nindex (kg/m2)", fc="#fbeae8", ec=RED, fs=8.4, lw=1.5)
arrow(axB, (0.13, 0.495), (0.455, 0.33), color=RED, lw=1.7)
arrow(axB, (0.545, 0.33), (0.87, 0.495), color=RED, lw=1.7)
axB.text(0.22, 0.375, "beta1 = +0.80\n(P = 3.2e-4)\nSteiger n.s.*", ha="center", va="center", fontsize=7.4, color=RED)
axB.text(0.66, 0.48, "beta2 = +0.361 / +0.593\n(P = 2e-32 / 8e-51)", ha="center", va="center", fontsize=7.6, color=RED)
axB.text(0.5, 0.225, "Indirect effect +0.291 / +0.477 (deleterious*)\nP = 5.8e-4 / 4.7e-4", ha="center", va="top", fontsize=8.2, color=RED, fontweight="bold")

# bottom-left MVMR note
box(axB, 0.02, 0.02, 0.60, 0.10,
    "Per-mediator MVMR (direct effects):\nadjusting SBP unmasks residual OR 1.44 / 1.90 (P = 0.002 / 5.6e-5)\nadjusting BMI moves to null: OR 0.89 / 1.06",
    fc="#fbfcfe", ec=NAVY, fs=7.6, lw=1.1)
axB.text(0.985, 0.015, "* BMI pathway fails Steiger\ndirectionality: instrument\nproperty, direction uncertain",
         ha="right", va="bottom", fontsize=6.8, color=GRAY, style="italic")

import os
os.makedirs("mr_results", exist_ok=True)
fig.savefig("mr_results/Figure1_triangulation_v2_300dpi.png", dpi=300, bbox_inches="tight")
fig.savefig("mr_results/Figure1_triangulation_v2_300dpi.tiff", dpi=300, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
print("saved")
