# -*- coding: utf-8 -*-
"""
plate_solver.io_export -- CSV / Excel export of solver results.

NEW capability (review-doc item): consolidated results used to be JSON-only
(the CHECKPOINT file). This adds a thin, dependency-light export layer on
top of that same summary structure so results can be opened directly in a
spreadsheet. It does not change what is computed or how -- purely a reporting
convenience, no SOLVER_VERSION concern.

Usage
-----
    from plate_solver.io_export import export_summary_csv, export_summary_xlsx
    export_summary_csv(summary, "results.csv")
    export_summary_xlsx(summary, "results.xlsx")   # requires openpyxl

`summary` is the same list-of-dict structure `compare_and_collect` /
`run_overnight`'s checkpoint already produce (or a list of plain dicts with
consistent keys) -- see validation.compare_and_collect for the exact shape.
"""
from __future__ import annotations
import csv
import json


def _flatten_rows(summary):
    """Normalize `summary` (list of dicts, possibly with list/tuple values)
    into a flat list of dicts with a stable, sorted column order."""
    rows = []
    for entry in summary:
        if isinstance(entry, dict):
            rows.append(entry)
        else:
            # Fallback: unknown row shape, keep it inspectable rather than
            # silently dropping data.
            rows.append({"value": repr(entry)})
    cols = []
    seen = set()
    for r in rows:
        for k in r.keys():
            if k not in seen:
                seen.add(k)
                cols.append(k)
    return rows, cols


def _stringify(v):
    if isinstance(v, (list, tuple, dict)):
        return json.dumps(v, default=str)
    return v


def export_summary_csv(summary, path):
    """Write `summary` (list of dicts) to a CSV file at `path`. Nested
    list/dict cell values are JSON-encoded into the cell (Excel/csv have no
    native nested-value support). Returns the number of rows written."""
    rows, cols = _flatten_rows(summary)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({k: _stringify(r.get(k, "")) for k in cols})
    return len(rows)


def export_summary_xlsx(summary, path):
    """Write `summary` to an .xlsx workbook at `path` using openpyxl.
    Raises ImportError with a clear message if openpyxl isn't installed
    (it is NOT one of this project's core dependencies -- mpmath/numpy/scipy
    are; this stays optional so a bare cluster environment doesn't need it
    just to run the solver)."""
    try:
        import openpyxl
    except ImportError as e:
        raise ImportError(
            "export_summary_xlsx requires openpyxl "
            "(pip install --break-system-packages openpyxl); "
            "use export_summary_csv for a dependency-free alternative."
        ) from e

    rows, cols = _flatten_rows(summary)
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.title = "results"
    ws.append(cols)
    for r in rows:
        ws.append([_stringify(r.get(k, "")) for k in cols])
    # Light formatting: bold header, freeze it, autosize-ish column widths.
    from openpyxl.styles import Font
    for cell in ws[1]:
        cell.font = Font(bold=True)
    ws.freeze_panes = "A2"
    for i, col in enumerate(cols, start=1):
        width = max(10, min(40, len(col) + 2))
        ws.column_dimensions[openpyxl.utils.get_column_letter(i)].width = width
    wb.save(path)
    return len(rows)


def load_checkpoint_summary(checkpoint_path):
    """Read a solver JSON checkpoint (CHECKPOINT env var target) and return
    its `summary` list, ready for export_summary_csv/xlsx. Returns [] if the
    file doesn't exist or has no `summary` key (rather than raising), since
    this is meant for ad-hoc post-run inspection."""
    try:
        with open(checkpoint_path, encoding="utf-8") as f:
            data = json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return []
    return data.get("summary", [])
