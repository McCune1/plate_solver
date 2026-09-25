# -*- coding: utf-8 -*-
"""
bundle_results_2026-09-24b.py -- collect every output of the combined
2026-09-24b queue into ONE zip for download (Python 3.6, stdlib only).

Run from Ansys/NewAnsys/ after the queue (the combined submit script calls it):
    python3 bundle_results_2026-09-24b.py [SLURM_JOB_ID]
Writes results_2026-09-24b_<jobid>.zip containing, for every deck in
ansys_queue_manifest_combined_2026-09-24b.txt: its _out.txt transcript and
every text output it wrote (modes / discriminator / ls* / freqs / samples ...),
plus the score_*.txt files, the targets JSONs, the manifest, the queue
.out/.err logs, and BUNDLE_INDEX.txt (what is in the zip, per deck, with
a MISSING line for anything expected but absent). Binary MAPDL files
(.rst, .db, .full, ...) are never included.
"""
import glob
import os
import re
import sys
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
MAN = "ansys_queue_manifest_combined_2026-09-24b.txt"
DATE_SFX = "_2026-09-24.inp"
EXTRA = ["score_p5mix_e31m_fe_2026-09-24.txt", "score_pcm_fe_2026-09-24.txt",
         "score_parashar_fe_2026-09-24.txt", "score_p5seg_e31m_fe_2026-09-24.txt",
         "targets_p5mix_e31m_2026-09-24.json", "targets_pcm_2026-09-24.json",
         "targets_parashar_2026-09-24.json", "targets_p5seg_e31m_2026-09-24.json",
         MAN]


CFOPEN = re.compile(r"^\*CFOPEN,([A-Za-z0-9_]+),(txt|csv)", re.M | re.I)


def deck_outputs(deck):
    """(expected, present): the deck's own *CFOPEN files plus its _out.txt."""
    src = open(os.path.join(HERE, deck), errors="replace").read()
    exp = [deck.replace(".inp", "_out.txt")]
    exp += ["%s.%s" % (n, e.lower()) for n, e in CFOPEN.findall(src)]
    exp = sorted(set(exp))
    return exp, [f for f in exp if os.path.exists(os.path.join(HERE, f))]


def main():
    jobid = sys.argv[1] if len(sys.argv) > 1 else "nojob"
    decks = [l.strip() for l in open(os.path.join(HERE, MAN))
             if l.strip() and not l.startswith("#")]
    index = ["results bundle for %s (job %s)" % (MAN, jobid), ""]
    members = set()
    for d in decks:
        exp, outs = deck_outputs(d)
        index.append("%s" % d)
        for o in exp:
            if o in outs:
                index.append("   %s" % o)
                members.add(o)
            else:
                index.append("   MISSING %s" % o)
    index.append("")
    for f in EXTRA + [os.path.basename(x) for x in
                      glob.glob(os.path.join(HERE, "ansys_queue_combined_2026-09-24b_*.out"))
                      + glob.glob(os.path.join(HERE, "ansys_queue_combined_2026-09-24b_*.err"))]:
        if os.path.exists(os.path.join(HERE, f)):
            members.add(f)
            index.append("extra: %s" % f)
        else:
            index.append("extra MISSING: %s" % f)
    zname = "results_2026-09-24b_%s.zip" % jobid
    with zipfile.ZipFile(os.path.join(HERE, zname), "w", zipfile.ZIP_DEFLATED) as z:
        for m in sorted(members):
            z.write(os.path.join(HERE, m), arcname="results_2026-09-24b/" + m)
        z.writestr("results_2026-09-24b/BUNDLE_INDEX.txt", "\n".join(index) + "\n")
    n_missing = sum(1 for l in index if "MISSING" in l)
    print("wrote %s: %d files, %d MISSING lines (see BUNDLE_INDEX.txt)"
          % (zname, len(members), n_missing))


if __name__ == "__main__":
    main()
