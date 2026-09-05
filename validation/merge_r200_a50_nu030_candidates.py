# -*- coding: utf-8 -*-
"""
merge_r200_a50_nu030_candidates.py -- no numerical work.
Combines r200_a50_nu030_capture_{oop,ip}.json into the shape
probe_ff_adjudicate_candidates.py expects. Submit via
submit_merge_r200_a50_nu030.sh (login node cannot run python).
ANSYS txt files are not required for the merge.
"""
import json
import sys

OOP_FILE = "r200_a50_nu030_capture_oop.json"
IP_FILE = "r200_a50_nu030_capture_ip.json"
OUT_FILE = "r200_a50_nu030_candidates.json"


def main():
    with open(OOP_FILE) as fh:
        oop = json.load(fh)
    with open(IP_FILE) as fh:
        ip = json.load(fh)

    if not oop.get("complete") or not ip.get("complete"):
        print("FATAL: one or both capture files are not marked complete. "
              f"OOP complete={oop.get('complete')}, "
              f"IP complete={ip.get('complete')}.")
        sys.exit(2)

    for key in ("r0_2b", "two_T_pi", "nu", "E", "rho"):
        if oop[key] != ip[key]:
            print(f"FATAL: geometry/material mismatch on {key}: "
                  f"OOP={oop[key]!r} IP={ip[key]!r}")
            sys.exit(2)

    merged = dict(
        tag="r200_a50_nu030",
        r0_2b=oop["r0_2b"], two_T_pi=oop["two_T_pi"], nu=oop["nu"],
        E=oop["E"], rho=oop["rho"],
        oop_n_dofs=oop["n_dofs"], oop_xmax=oop["xmax"],
        ip_n_dofs=ip["n_dofs"], ip_xmax=ip["xmax"],
        fe_oop_file="geomsweep_oop_r200_a50_nu030.txt",
        fe_ip_file="geomsweep_ip_r200_a50_nu030.txt",
        fe_oop_rigid_count=6, fe_ip_rigid_count=3,
        source=f"merged from {OOP_FILE} (job {oop['job_id']}) + "
               f"{IP_FILE} (job {ip['job_id']}), rank-16 crossed "
               "geometry r0/2b=2.0 2T/pi=0.5 nu=0.30",
        oop_omegas=oop["omegas"],
        ip_omegas=ip["omegas"],
    )
    with open(OUT_FILE, "w") as fh:
        json.dump(merged, fh, indent=1)
    print(f"wrote {OUT_FILE}: {len(merged['oop_omegas'])} OOP + "
          f"{len(merged['ip_omegas'])} IP candidates.")
    print("Next: sbatch submit_ff_adjudicate_r200_a50_nu030.sh")
    print("(only after Ansys/NewAnsys/geomsweep_{oop,ip}_r200_a50_nu030.txt "
          "exist — they already do).")


if __name__ == "__main__":
    main()
