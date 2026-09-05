# -*- coding: utf-8 -*-
"""Rebuild PLATE_SOLVER_REPO_DIGEST.md from github_repo/.

The digest is the single-file review bundle at the project root (one
directory above github_repo/). Binary files are listed, not inlined.
ANSYS MAPDL *_out.txt logs are condensed to their frequency tables.

Run from anywhere:
    python github_repo/scripts/build_repo_digest.py
"""
from __future__ import annotations

import datetime as _dt
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]          # github_repo/
PROJECT = ROOT.parent                               # Plate_Solver_Package/
OUT = PROJECT / "PLATE_SOLVER_REPO_DIGEST.md"

SKIP_DIR_NAMES = {"__pycache__", ".git", ".ruff_cache"}
SKIP_SUFFIXES_BINARY = {".pdf", ".png", ".jpg", ".jpeg", ".gif"}
SKIP_SUFFIXES_DEBRIS = {".pyc", ".aux", ".fls", ".fdb_latexmk", ".synctex.gz"}
SKIP_LATEX_DEBRIS_NAMES = {
    "PAPER1_FREEFREE_DRAFT.log",
    "PAPER1_FREEFREE_DRAFT.out",
    "PAPER1_FREEFREE_DRAFT.aux",
    "PAPER1_FREEFREE_SUPPLEMENTARY.log",
    "PAPER1_FREEFREE_SUPPLEMENTARY.out",
    "PAPER1_FREEFREE_SUPPLEMENTARY.aux",
}

FREQ_START = re.compile(
    r"(FREQUENCIES FROM BLOCK|FREQUENCY \(HERTZ\)|"
    r"\*\*\* FREQUENCIES|\bMODE\s+FREQUENCY)",
    re.I,
)
FREQ_LINE = re.compile(
    r"^\s*\d+\s+[-\d.Ee+]+",
)


def rel(p: Path) -> str:
    return p.relative_to(ROOT).as_posix()


def iter_repo_files():
    for dirpath, dirnames, filenames in os.walk(ROOT):
        dirnames[:] = [d for d in sorted(dirnames) if d not in SKIP_DIR_NAMES]
        for name in sorted(filenames):
            p = Path(dirpath) / name
            if name in SKIP_LATEX_DEBRIS_NAMES:
                continue
            if p.suffix.lower() in SKIP_SUFFIXES_DEBRIS:
                continue
            yield p


def is_binary_path(p: Path) -> bool:
    return p.suffix.lower() in SKIP_SUFFIXES_BINARY


def kb(n: int) -> str:
    return f"{max(1, round(n / 1024))} KB"


def condense_ansys_log(text: str, path: str, nbytes: int) -> str:
    lines = text.splitlines()
    blocks: list[str] = []
    i = 0
    while i < len(lines):
        if FREQ_START.search(lines[i]):
            chunk = [lines[i]]
            i += 1
            blanks = 0
            while i < len(lines):
                line = lines[i]
                if FREQ_START.search(line) and len(chunk) > 3:
                    break
                chunk.append(line)
                if line.strip() == "":
                    blanks += 1
                    if blanks >= 3 and any(FREQ_LINE.match(x) for x in chunk):
                        i += 1
                        break
                else:
                    blanks = 0
                # stop if we wander into unrelated MAPDL dump
                if len(chunk) > 400:
                    break
                i += 1
            blocks.append("\n".join(chunk).rstrip())
        else:
            i += 1
    if not blocks:
        return (
            f"[Condensed: {path} is a raw ANSYS MAPDL solver log "
            f"({kb(nbytes)}). No frequency table could be extracted "
            f"automatically; the full raw log is in the repository.]\n"
        )
    header = (
        f"[Condensed: this file is a raw ANSYS MAPDL solver log "
        f"({kb(nbytes)}). License banner, element/node listings, and "
        f"solver performance statistics are omitted here as boilerplate "
        f"with no review-relevant content; the full raw log is at "
        f"`{path}` in the repository. Modal frequency table(s) below "
        f"are extracted verbatim.]\n"
    )
    return header + "\n\n".join(blocks) + "\n"


def should_condense_ansys(p: Path) -> bool:
    if "ansys" not in rel(p).split("/"):
        return False
    name = p.name.lower()
    return name.endswith("_out.txt") or name.endswith("_out.log")


def read_text(p: Path) -> str:
    raw = p.read_bytes()
    for enc in ("utf-8", "utf-8-sig", "cp1252", "latin-1"):
        try:
            return raw.decode(enc)
        except UnicodeDecodeError:
            continue
    return raw.decode("utf-8", errors="replace")


def file_body(p: Path) -> str:
    r = rel(p)
    n = p.stat().st_size
    text = read_text(p)
    if should_condense_ansys(p):
        return condense_ansys_log(text, r, n)
    if not text.endswith("\n"):
        text += "\n"
    return text.replace("\r\n", "\n")


def build() -> str:
    files = list(iter_repo_files())
    text_files = [p for p in files if not is_binary_path(p)]
    bin_files = [p for p in files if is_binary_path(p)]
    today = _dt.date.today().isoformat()

    intro = f"""# plate_solver — Consolidated Repository Reference

**Generated {today} from `github_repo/`.** Live Paper 1 manuscript:
`github_repo/paper/PAPER1_FREEFREE_DRAFT.tex` +
`github_repo/paper/PAPER1_FREEFREE_SUPPLEMENTARY.tex` (compiled PDFs next
to those sources). See `Project Knowledge/CURRENT_PAPER.md`. This digest
is a review bundle; if it disagrees with `github_repo/paper/`, the paper
directory wins. Rebuild with
`python github_repo/scripts/build_repo_digest.py`.

This is a single-file, text-only bundle of the `plate_solver` repository
(the companion code/data for "Exact Wave-Function Formulation and
Calibrated Discrimination Instruments for Completely Free Annular Sector
Plates", G. McCune and D. Stutts, Missouri University of Science and
Technology, draft). It exists so the repository and the paper draft can
be reviewed together by pasting one document into a chat, without needing
repository access.

It contains, in full, every text file in the repository: the `plate_solver`
package source, its test suite, the validation/reproduction scripts and
their evidence logs, the SLURM launch scripts, the ANSYS APDL decks, the
checkpoint data, and the paper's LaTeX source (main text and Supplementary
Material). The binary file types in the repository (compiled PDF, PNG
figures) cannot be embedded in a text file and are listed below, along
with their repository paths. Several large raw ANSYS solver logs (license
banner, element/node listings, solver performance statistics — all
boilerplate) are condensed to just their modal-frequency result tables;
each condensation says so and gives the path to the full original.

**How this file is organized:** below the file tree, every text file's full
contents follow, each preceded by a `FILE: <relative path>` divider line —
paths are relative to the repository root and match the tree exactly.
Markdown files (README.md, etc.) appear as their raw markdown; code/config
files appear as raw source. There is intentionally no code-fencing wrapper
around each entry (several of the source READMEs contain their own fenced
code examples, which would otherwise break nesting) — the `FILE:` divider
lines are the section boundaries.

**For a reviewer using this file:** the paper draft is in the file tree at
`paper/PAPER1_FREEFREE_DRAFT.tex`, with its Supplementary Material at
`paper/PAPER1_FREEFREE_SUPPLEMENTARY.tex` and a standalone graphical
abstract at `paper/GRAPHICAL_ABSTRACT.pdf` (binary, listed not inlined).
The top-level `README.md` (also below) explains the repository's structure
and status. `plate_solver/` is the package the paper's results were
computed with; `validation/` and `data/checkpoints/` are what backs the
paper's validated tables and figures, including `validation/geomsweep/`
and `ansys/geomsweep/` (the 48-geometry sweep, current headline 1256/1475
within 3%) and `validation/ip_mac/` (the §6.6 MAC bar). `SOLVER_VERSION`
in `plate_solver/config.py` is `2026-07-10.s10`.

## File tree

```
"""
    tree_lines = [rel(p) for p in files]
    intro += "\n".join(tree_lines) + "\n```\n\n"
    intro += "## Binary files not included (see repository for these)\n\n"
    for p in bin_files:
        intro += f"- `{rel(p)}`  ({kb(p.stat().st_size)})\n"
    intro += "\n---\n\n# File contents\n\n"

    parts = [intro]
    for p in text_files:
        r = rel(p)
        body = file_body(p)
        parts.append(
            "==============================================================================\n"
            f"FILE: {r}\n"
            "==============================================================================\n\n"
            f"{body}"
        )
        if not body.endswith("\n"):
            parts[-1] += "\n"
        parts.append("\n")
    return "".join(parts)


def main() -> None:
    text = build()
    OUT.write_text(text, encoding="utf-8", newline="\n")
    n_file = text.count("\nFILE: ")
    print(f"wrote {OUT}  ({len(text):,} bytes, {n_file} FILE: sections)")


if __name__ == "__main__":
    main()
