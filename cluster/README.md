# cluster/

SLURM launch scripts for running the full `plate_solver` batch as it was
actually run for this paper (Missouri S&T's "Mill" HPC cluster, AMD EPYC
7502 32-core nodes). These are cluster-specific — the `#SBATCH` header
(partition name, core count, wall-clock limit, module-load lines) will need
to be adapted to whatever scheduler and Python environment you're using.
Every environment variable / flag they set is documented in
`plate_solver/cli.py`'s argparse layer.

| File | Runs |
|---|---|
| `submit_plate_solver.sh` | the full overnight batch: cantilever regression spot-check, free-free spectra sweep (both radial edges free, parallel across cores), literature/FE comparison, orthotropic-material sanity check |
| `submit_rect_ip.sh` | rectangular in-plane cantilever validation |
| `submit_diag.sh` | diagnostic/short-run mode for checking a deployment before committing to a full batch |

For a single-point sanity check without a scheduler, run
`python3 -m plate_solver.cli` directly (see the top-level README's
Quickstart) — but note that a real mode search does not finish quickly on a
laptop; these scripts exist because the full sweep needs the parallelism.

Each script also hardcodes a venv activation path from the original cluster
deployment (`/home/ghmkfh/PythonMill/Plate_Solver_Package/venv`) — point
this at your own environment before submitting.
