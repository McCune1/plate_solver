"""Geometry schematic and identifiability-vs-thickness figure for Paper 4.
2026-09-24: values at e31 = -4.1 (LESSONS Sec 18.246)."""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, FancyArrowPatch
import numpy as np

OUT = Path("paper_figures")

# --- geometry schematic ---
fig, ax = plt.subplots(figsize=(6.6, 3.4))
ax.set_xlim(0, 10)
ax.set_ylim(0, 5.2)
ax.set_aspect("equal")
ax.axis("off")

# radial slice: host and two piezo layers
# x: r from ri to ro, y: z
x0, w = 1.2, 6.5
h_host, h_p = 1.35, 0.55
y_host = 1.7
ax.add_patch(Rectangle((x0, y_host), w, h_host, facecolor="#d9d9d9",
                       edgecolor="k", lw=1.0, label="host"))
ax.add_patch(Rectangle((x0, y_host + h_host), w, h_p, facecolor="#6baed6",
                       edgecolor="k", lw=1.0))
ax.add_patch(Rectangle((x0, y_host - h_p), w, h_p, facecolor="#6baed6",
                       edgecolor="k", lw=1.0))
# electrodes: inner (interface) solid, outer dashed for OC
for y in (y_host, y_host + h_host):
    ax.plot([x0, x0 + w], [y, y], color="#222", lw=2.2)
for y in (y_host - h_p, y_host + h_host + h_p):
    ax.plot([x0, x0 + w], [y, y], color="#c0392b", lw=2.2, ls="--")

ax.annotate("", xy=(x0, 0.55), xytext=(x0 + w, 0.55),
            arrowprops=dict(arrowstyle="<->", color="k", lw=0.9))
ax.text(x0 + w / 2, 0.22, r"$r_i\longrightarrow r_o$", ha="center", fontsize=10)
ax.annotate("", xy=(0.75, y_host), xytext=(0.75, y_host + h_host),
            arrowprops=dict(arrowstyle="<->", color="k", lw=0.9))
ax.text(0.35, y_host + h_host / 2, r"$2h$", ha="center", va="center",
        fontsize=10, rotation=90)
ax.annotate("", xy=(8.05, y_host + h_host), xytext=(8.05, y_host + h_host + h_p),
            arrowprops=dict(arrowstyle="<->", color="k", lw=0.9))
ax.text(8.45, y_host + h_host + h_p / 2, r"$h_1$", ha="left", va="center",
        fontsize=10)

ax.text(x0 + w / 2, y_host + h_host / 2, "steel host", ha="center", va="center",
        fontsize=10)
ax.text(x0 + w / 2, y_host + h_host + h_p / 2, "PZT", ha="center", va="center",
        fontsize=9, color="k")
ax.text(x0 + w / 2, y_host - h_p / 2, "PZT", ha="center", va="center",
        fontsize=9, color="k")

ax.plot([9.1, 9.7], [4.55, 4.55], color="#222", lw=2.2)
ax.text(9.85, 4.55, "inner electrodes (grounded)", va="center", fontsize=8)
ax.plot([9.1, 9.7], [4.15, 4.15], color="#c0392b", lw=2.2, ls="--")
ax.text(9.85, 4.15, "outer electrodes (bus, $Q=0$ if OC)", va="center", fontsize=8)

ax.set_xlim(0, 13.2)
fig.tight_layout()
fig.savefig(OUT / "piezo_p4_geometry.png", dpi=160, bbox_inches="tight")
plt.close()

# --- identifiability vs thickness ---
ratio = np.array([1 / 12, 1 / 8, 1 / 5])
unc_sc = np.array([330.26, 109.21, 32.69])  # e31=-4.1 cc_sc_e31_e31m (+4.1: 888.5, 293.8, 87.9)  # 2026-09-23 Sec 18.231 (duan: 730.8, 241.6, 72.3)
unc_oc = np.array([0.8701, 0.6280, 0.4526])  # e31=-4.1 h1_ratio_sweep_e31m G2 (+4.1: 2.133, 1.479, 1.004)
k2 = np.array([0.037881, 0.053296, 0.075719])  # e31=-4.1 (+4.1: 0.005573, 0.008054, 0.011912)

fig, ax = plt.subplots(figsize=(6.4, 3.8))
ax.semilogy(ratio, unc_sc, "s-", color="#2166ac", lw=1.4, ms=7,
            label=r"C–C short-circuit $e_{31}$")
ax.semilogy(ratio, unc_oc, "o-", color="#b2182b", lw=1.4, ms=7,
            label=r"F–F open-circuit $e_{31}$")
ax.set_xticks(ratio)
ax.set_xticklabels([r"$1/12$", r"$1/8$", r"$1/5$"])
ax.set_xlabel(r"$h_1/2h$")
ax.set_ylabel("rel. uncertainty in $e_{31}$ at 0.01% freq.")
ax.set_ylim(0.2, 1000)
ax.legend(frameon=False, fontsize=9)
ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.7)
fig.tight_layout()
fig.savefig(OUT / "piezo_p4_identifiability.png", dpi=160)
plt.close()
print("wrote geometry + identifiability")
