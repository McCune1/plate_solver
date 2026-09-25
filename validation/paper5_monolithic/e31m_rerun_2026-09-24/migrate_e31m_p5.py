# Writes *_e31m_2026-09-24.py copies of the Paper 5 probes at e31 = -4.1.
import re, os
BANNER = ("# e31 = -4.1 COPY (LESSONS Sec 18.246) of {parent}.\n"
          "# Written by migrate_e31m_p5.py. Only e31 changes, plus the windows that\n"
          "# bracket the F-F SC root (471.36 -> 519.84 rad/s at -4.1, i.e. every\n"
          "# window near that root shifts by +48.47), the e31 inverse bracket\n"
          "# [1,8] -> [-8,-1] (wrong-sign control [-8,-1] -> [1,8]; the twin\n"
          "# 2*8.95 - e31 moves to 22.0), relative errors divide by |E31_TRUE|,\n"
          "# and output files carry the suffix _e31m. Fixed before running.\n")
EDITS = {
  "probe_piezo_p5_cf_boundary_2026-09-20.py": [
      ("E31, E33 = 4.1, 14.1", "E31, E33 = -4.1, 14.1")],
  "probe_piezo_p5_e31_inverse_2026-09-20.py": [
      ("E31_TRUE, E33 = 4.1, 14.1", "E31_TRUE, E33 = -4.1, 14.1"),
      ("E31_LO, E31_HI = 1.0, 8.0", "E31_LO, E31_HI = -8.0, -1.0"),
      ("NEG_LO, NEG_HI = -8.0, -1.0", "NEG_LO, NEG_HI = 1.0, 8.0"),
      ("g0_rel = abs(e31_hat - E31_TRUE) / E31_TRUE", "g0_rel = abs(e31_hat - E31_TRUE) / abs(E31_TRUE)")],
  "probe_piezo_p5_n0_overtones_2026-09-20.py": [
      ("E31_TRUE, E33 = 4.1, 14.1", "E31_TRUE, E33 = -4.1, 14.1")],
  "probe_piezo_p5_thickness_sweep_2026-09-20.py": [
      ("E31_TRUE, E33 = 4.1, 14.1", "E31_TRUE, E33 = -4.1, 14.1")],
  "probe_piezo_p5_yomega_sense_2026-09-20.py": [
      ("E31, E33 = 4.1, 14.1", "E31, E33 = -4.1, 14.1"),
      ("OMEGA_NEAR_POLE = 473.0", "OMEGA_NEAR_POLE = 521.5"),
      ("CP_FF_LO, CP_FF_HI = 460.0, 490.0", "CP_FF_LO, CP_FF_HI = 505.0, 535.0"),
      ("for omega_s in [460.0, 468.0, 471.0, omega_coupled - 0.5,",
       "for omega_s in [508.5, 516.5, 519.5, omega_coupled - 0.5,"),
      ("omega_coupled + 0.5, 478.0, 485.0]:", "omega_coupled + 0.5, 526.5, 533.5]:")],
}
for f, eds in EDITS.items():
    s = open(f).read()
    for a, b in eds:
        assert s.count(a) >= 1, (f, a)
        s = s.replace(a, b)
    s, k = re.subn(r'(_results)(\.json)', r'\1_e31m\2', s)
    out = f.replace("2026-09-20.py", "e31m_2026-09-24.py")
    open(out, "w", newline="\n").write(BANNER.format(parent=f) + s)
    print(out, "json", k)
