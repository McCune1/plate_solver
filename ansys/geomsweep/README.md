# ansys/geomsweep/

Independent finite-element benchmark for `paper/PAPER1_FREEFREE_DRAFT.tex`
§6.4 — a second, isolated material set (E=210 GPa, nu=0.35, rho=7800
kg/m^3, H=0.08 m; deliberately not the nu=0.30 set behind the main §6/§6.3
benchmark) run at 4 radius ratios (r0/2b = 1.5, 1.66667, 2.0, 2.5) across
all 6 sector angles each, both motion types. Binary ANSYS scratch/result
files are excluded, as in `ansys/` — see its README and `.gitignore`.

| Deck | Used for |
|---|---|
| `ansys_geomsweep_oop_r150.inp` / `_r167.inp` / `_r200.inp` / `_r250.inp` | flexural (SHELL281) FE benchmark at each radius ratio, all 6 angles per deck, §6.4 |
| `ansys_geomsweep_ip_r150.inp` / `_r167.inp` / `_r200.inp` / `_r250.inp` | extensional (PLANE183) FE benchmark at each radius ratio, all 6 angles per deck, §6.4 |
| `generate_geomsweep_ansys_decks.py` | generator for all 8 decks above — do not hand-edit the `.inp` files, edit this and regenerate instead |
| `submit_ansys_geomsweep_oop_r*.sh` / `_ip_r*.sh` | SLURM wrappers for the 8 decks; adapt the SBATCH header, ANSYS module, and hardcoded working directory to your own scheduler and install location |

Each `.inp` has a matching `*_out.txt`, the raw SLURM-captured MAPDL solver
log for all 6 angles in that deck (license banner, element/node listings,
solver performance statistics included — this is the unprocessed job
output, kept for provenance). **The actual clean per-angle results are in
the separate `geomsweep_{oop,ip}_r*_a*.txt` files** (48 total, 6 angles ×
4 ratios × 2 parts) — each deck opens its own output file via `*CFOPEN`
and `*VWRITE`s a `mode  f[Hz]  Omega_lit` table into it directly, so this
data does *not* appear in the `_out.txt` solver log at all (verified: the
Omega_lit values are absent from the raw logs, present only in these small
per-angle files). These are the actual FE ground truth §6.4's matching
analysis used.

**Reproducibility note (CLOSED 2026-07-23, cluster-confirmed):** the
original script that performed the python-vs-FE nearest-neighbor matching
(producing §6.4's original "1232/1466 candidates (84.0%) matched within
3%" claim) was run interactively on 2026-07-22 and never saved — an
external review flagged this as a real reproducibility gap.
`validation/geomsweep/probe_geomsweep_fe_match.py` is a from-scratch
reconstruction from these FE files plus the python-side manifests, run
for real on the cluster as job 2328808
(`validation/geomsweep/geomsweep_fe_match_2328808.out`): n_total 1474,
matched 1231/83.5% (mean 0.577%, worst 3.000%) against the original
1466/84.0% (mean 0.58%, worst 3.00%). That recapture is itself superseded:
the live paper (§6.4 and SM Table S.3) reports 1256/1475 (85.2%, mean
0.574%, worst 2.990%), with a footnote listing 1231/1474 as an intermediate
figure. Job 2416746's match log and CSV are in
`validation/geomsweep/geomsweep_fe_match_s10v2_2416746.out`. Do not treat
1231/1474 as the manuscript's current headline. The extensional-family
contamination tolerance that could not be recovered from the original
interactive run is documented in `validation/geomsweep/README.md`.
Rigid-body modes are handled per §6.3's stated procedure.
