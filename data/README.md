# data/checkpoints/

Checkpoint JSON files written by `plate_solver.cli` / `run_overnight.py`
during the overnight cantilever-regression + free-free-sweep batch. These
are the inputs the `scripts/probe_regen_*_figs.py` figure scripts read.

| File | Corresponds to |
|---|---|
| `research23_checkpoint.json` | main run |
| `research23_checkpoint_p1tail.json` | Part-1 tail continuation |
| `research23_checkpoint_p2.json` | Part-2 continuation |

Each checkpoint records the `SOLVER_VERSION` it was produced under
(`plate_solver/config.py`); a checkpoint is only meaningful for the solver
version it was written with — do not mix a checkpoint from one version with
a script run against another.
