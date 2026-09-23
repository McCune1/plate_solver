"""Figures for Paper 5 (monolithic piezoelectric ring)."""
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
import numpy as np

OUT = Path("paper_figures")
OUT.mkdir(exist_ok=True)

# --- 1. geometry schematic: single homogeneous ceramic, full thickness 2H ---
fig, ax = plt.subplots(figsize=(6.6, 3.2))
ax.set_xlim(0, 10.6)
ax.set_ylim(0, 5.0)
ax.set_aspect("equal")
ax.axis("off")

x0, w = 1.2, 6.5
h_ceramic = 1.8
y0 = 1.6
ax.add_patch(Rectangle((x0, y0), w, h_ceramic, facecolor="#6baed6",
                        edgecolor="k", lw=1.0))
ax.text(x0 + w / 2, y0 + h_ceramic / 2, "homogeneous PZT-4\n(thickness-poled)",
        ha="center", va="center", fontsize=9.5)

for y in (y0, y0 + h_ceramic):
    ax.plot([x0, x0 + w], [y, y], color="#222", lw=2.4)

ax.annotate("", xy=(x0 + w * 0.15, y0 + h_ceramic + 0.35),
            xytext=(x0 + w * 0.15, y0 - 0.35),
            arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=1.6))
ax.text(x0 + w * 0.15 - 0.55, y0 + h_ceramic / 2, "$P$", ha="center",
        va="center", fontsize=10, color="#c0392b")

ax.annotate("", xy=(x0, 0.55), xytext=(x0 + w, 0.55),
            arrowprops=dict(arrowstyle="<->", color="k", lw=0.9))
ax.text(x0 + w / 2, 0.22, r"$r_i\longrightarrow r_o$", ha="center", fontsize=10)

ax.annotate("", xy=(0.75, y0), xytext=(0.75, y0 + h_ceramic),
            arrowprops=dict(arrowstyle="<->", color="k", lw=0.9))
ax.text(0.35, y0 + h_ceramic / 2, r"$2H$", ha="center", va="center",
        fontsize=10, rotation=90)

ax.plot([x0, x0 + 0.6], [4.55, 4.55], color="#222", lw=2.4)
ax.text(x0 + 0.8, 4.55, "face electrodes: grounded (SC) or floating bus, $Q=0$ (OC)",
        va="center", fontsize=8)

# Render at a generous fixed size (no 'tight'), then post-crop to the drawn
# content with an equal pixel pad on every side. matplotlib's bbox_inches="tight"
# was including the invisible Axes background box (sized by the wide xlim) rather
# than just the drawn artists, which padded the right side of this figure far more
# than the left; this post-crop sidesteps that entirely and is exactly symmetric.
import io
from PIL import Image
buf = io.BytesIO()
fig.savefig(buf, format="png", dpi=160)
plt.close()
buf.seek(0)
img = Image.open(buf).convert("RGB")
arr = np.array(img)
nonwhite = np.any(arr < 250, axis=2)
cols = np.where(nonwhite.any(axis=0))[0]
rows = np.where(nonwhite.any(axis=1))[0]
pad = 14
left, right = max(cols.min() - pad, 0), min(cols.max() + pad, arr.shape[1] - 1)
top, bottom = max(rows.min() - pad, 0), min(rows.max() + pad, arr.shape[0] - 1)
img.crop((left, top, right + 1, bottom + 1)).save(OUT / "piezo_p5_geometry.png")

# --- 2. Y_sense resonance sweep (9-point pre-registered gate sweep, no new compute) ---
# 2026-09-23 (LESSONS Sec 18.231): consistent-projection rerun of the
# pre-registered sweep, piezo_p5_yomega_sense_results.json.
sweep = [
    (460.0, -7.673041357247039e-11),
    (468.0, -2.809962953233649e-10),
    (470.8646021920983, -1.9436000771931203e-09),
    (471.0, -2.668790293105272e-09),
    (471.3546021920983, -9.763123604459714e-08),
    (471.3746021920983, 9.764965584198399e-08),
    (471.8646021920983, 1.9620198745501872e-09),
    (478.0, 1.5637479197710927e-10),
    (485.0, 8.084682698928616e-11),
]
om = np.array([p[0] for p in sweep])
Y = np.array([p[1] for p in sweep])
pole = 471.3646021920983

fig, axes = plt.subplots(1, 2, figsize=(9.4, 3.6))

ax = axes[0]
ax.semilogy(om, np.abs(Y), "o-", color="#2166ac", lw=1.3, ms=5)
ax.axvline(pole, color="#c0392b", ls="--", lw=1.0)
ax.set_xlabel(r"$\omega\ (\mathrm{rad\,s^{-1}})$")
ax.set_ylabel(r"$|Y_{\mathrm{sense}}[r_*]|$")
ax.set_title("(a) full pre-registered sweep")
ax.grid(True, which="both", ls=":", lw=0.5, alpha=0.6)

ax = axes[1]
near = om[np.abs(om - pole) < 1.0]
Yn = Y[np.abs(om - pole) < 1.0]
colors = ["#2166ac" if v < 0 else "#b2182b" for v in Yn]
ax.scatter(near, Yn, c=colors, s=40, zorder=3)
ax.plot(near, Yn, color="#888", lw=1.0, zorder=2)
ax.axhline(0, color="k", lw=0.6)
ax.axvline(pole, color="#c0392b", ls="--", lw=1.0, label=r"$\omega=471.3646$")
ax.set_xlabel(r"$\omega\ (\mathrm{rad\,s^{-1}})$")
ax.set_ylabel(r"$Y_{\mathrm{sense}}[r_*]$")
ax.set_title("(b) sign flip through the pole")
ax.legend(frameon=False, fontsize=8, loc="upper left")
ax.grid(True, ls=":", lw=0.5, alpha=0.6)

fig.tight_layout()
fig.savefig(OUT / "piezo_p5_ysense_sweep.png", dpi=160)
plt.close()

# --- 3. stiffening vs mechanical BC ---
bcs = ["F-F", "C-C", "C-F", "F-C"]
splits = [1.939229678534127, 1.783814410565254, 2.023078154208058, 1.439419973002432]  # Sec 18.231
colors = ["#4393c3", "#4393c3", "#92c5de", "#92c5de"]

fig, ax = plt.subplots(figsize=(5.2, 3.4))
bars = ax.bar(bcs, splits, color=colors, edgecolor="k", lw=0.8, width=0.55)
for b, v in zip(bars, splits):
    ax.text(b.get_x() + b.get_width() / 2, v + 0.03, f"{v:.2f}%",
            ha="center", va="bottom", fontsize=9)
ax.set_ylabel(r"SC-vs-elastic split, $n=0$ fundamental (\%)")
ax.set_ylim(0, 2.9)
ax.grid(True, axis="y", ls=":", lw=0.5, alpha=0.6)
fig.tight_layout()
fig.savefig(OUT / "piezo_p5_bc_stiffening.png", dpi=160)
plt.close()

print("wrote piezo_p5_geometry.png, piezo_p5_ysense_sweep.png, piezo_p5_bc_stiffening.png")
