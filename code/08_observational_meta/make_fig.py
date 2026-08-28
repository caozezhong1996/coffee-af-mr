# -*- coding: utf-8 -*-
"""Forest plot for the observational dose-response wing (Figure 4 candidate)."""
import sys, json
from pathlib import Path
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("percup_trends.csv")
res = json.load(open("meta_wing_results.json"))["percup_RE"]

fig, ax = plt.subplots(figsize=(7.2, 3.4))
ys = list(range(len(df), 0, -1))
for (_, r), y in zip(df.iterrows(), ys):
    ax.plot([r["lo"], r["hi"]], [y, y], "-", color="#444", lw=1.4)
    ax.plot(r["HR"], y, "s", color="#1f6fb2", ms=7)
    ax.text(0.915, y, r["study"], ha="left", va="center", fontsize=9)
    ax.text(1.10, y, f"{r['HR']:.3f} ({r['lo']:.3f}–{r['hi']:.3f})", ha="left", va="center", fontsize=9)
# pooled diamond
yp = 0
c = res
ax.fill([c["lo"], c["HR_per_cup_RE"], c["hi"], c["HR_per_cup_RE"]], [yp, yp+0.32, yp, yp-0.32],
        color="#c0392b", alpha=0.85)
ax.text(0.915, yp, "Random-effects pooled", ha="left", va="center", fontsize=9, fontweight="bold")
ax.text(1.10, yp, f"{c['HR_per_cup_RE']:.3f} ({c['lo']:.3f}–{c['hi']:.3f})", ha="left", va="center",
        fontsize=9, fontweight="bold")
ax.axvline(1.0, color="#888", ls="--", lw=0.8)
ax.set_xlim(0.93, 1.13)
ax.set_ylim(-0.8, len(df)+0.8)
ax.set_yticks([])
ax.set_xscale("log")
from matplotlib.ticker import NullLocator, FixedLocator, FixedFormatter
ax.xaxis.set_minor_locator(NullLocator())
ax.xaxis.set_major_locator(FixedLocator([0.95, 1.0, 1.05, 1.1]))
ax.xaxis.set_major_formatter(FixedFormatter(["0.95", "1.00", "1.05", "1.10"]))
ax.set_xlabel("Hazard ratio of incident AF per additional cup/day (log scale)")
ax.set_title(f"Observational dose-response wing: per-cup trend, 5 prospective cohorts "
             f"(I$^2$ = {res['I2_pct']:.0f}%)", fontsize=10)
fig.savefig("forest_meta_wing.png", bbox_inches="tight", dpi=300)
print("saved forest_meta_wing.png")
