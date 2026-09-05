# -*- coding: utf-8 -*-
"""
plate_solver -- Annular & Rectangular Sector Plate Vibration Solver.

Modularized 2026-07-01 from Research50.py (annular, orthotropy-validated)
and rect_int.py (rectangular). See LESSONS_LEARNED.md Sec. 22 for the full
provenance/rationale and what has/has not been re-verified post-split.

Import order matters: config MUST be imported first (it single-threads BLAS
before numpy/scipy ever get imported anywhere in the process).
"""
from __future__ import annotations

from .config import CONFIG, SolverConfig, SOLVER_VERSION

from .geometry import (
    PlateGeometry, MaterialModel, IsotropicMaterial, OrthotropicMaterial,
    FGMPlateProperties, RadialFGMMaterial, make_geometry,
)
from .boundary import (
    EdgeKind, EdgeSpec, BoundaryCondition, ClampedFreeOOP, ClampedFreeIP,
    FreeFreeOOP, FreeFreeIP, make_bc,
)
from .dispersion import (
    ExactEdgeSolver, cutoff_frequencies_part1, cutoff_frequencies_part2,
    RectangularCartesianOOP, RectangularCartesianIP,
)
from .detectors import (
    full_search, track, select_fill, find_modes_sigmin, sigma_min_from_K,
    equilibrated_logdet, equilibrated_nullvec_mp,
)
from .core_solvers import (
    OutOfPlaneSolver, InPlaneSolver, RectOOPAssembler, RectIPAssembler,
    AnnulusRadialOOP, AnnulusRadialIP,
)
from .plotting import plot_dispersion_curves, plot_mode_shape
from .validation import (
    compare_and_collect, print_summary, run_shi_validation,
    run_mcgee_validation, run_rect_validation, run_rect_ip_validation,
    PAPER_PART1, PAPER_PART2,
)
from .run_overnight import run_overnight

__all__ = [
    "CONFIG", "SolverConfig", "SOLVER_VERSION",
    "PlateGeometry", "MaterialModel", "IsotropicMaterial", "OrthotropicMaterial",
    "FGMPlateProperties", "RadialFGMMaterial", "make_geometry",
    "EdgeKind", "EdgeSpec", "BoundaryCondition", "ClampedFreeOOP", "ClampedFreeIP",
    "FreeFreeOOP", "FreeFreeIP", "make_bc",
    "ExactEdgeSolver", "cutoff_frequencies_part1", "cutoff_frequencies_part2",
    "RectangularCartesianOOP", "RectangularCartesianIP",
    "full_search", "track", "select_fill", "find_modes_sigmin", "sigma_min_from_K",
    "equilibrated_logdet", "equilibrated_nullvec_mp",
    "OutOfPlaneSolver", "InPlaneSolver", "RectOOPAssembler", "RectIPAssembler",
    "AnnulusRadialOOP", "AnnulusRadialIP",
    "plot_dispersion_curves", "plot_mode_shape",
    "compare_and_collect", "print_summary", "run_shi_validation",
    "run_mcgee_validation", "run_rect_validation", "run_rect_ip_validation",
    "PAPER_PART1", "PAPER_PART2", "run_overnight",
]
