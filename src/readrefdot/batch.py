"""readrefdot-batch - run readrefdot over the rows of a TSV manifest.

Each row names its own inputs, output location and per-plot arguments, so one file
can mix samples, references and settings. Rows whose output already exists are
skipped, making the wrapper safe to re-run after adding rows or after a failure.

Required columns  bam, readid, reference, outdir, suffix
Optional columns  k, min-seg, merge-gap, panel-mm, colour_main, colour_ext,
                  ref-lines, read-lines, monomer, monomer-period, monomer-cut
                  (blank means "use the default"; '-' spellings also accepted with '_')

`suffix` is the output file stem: a row writes <outdir>/<suffix>.quad.png and .pdf.
"""

import argparse
import csv
import os
import sys
from collections import OrderedDict

from .annotate import Lines
from .plot import Params, quad
from .read import ReadNotFound, load

REQUIRED = ("bam", "readid", "reference", "outdir", "suffix")
OPTIONAL = ("k", "min-seg", "merge-gap", "panel-mm", "colour_main", "colour_ext",
            "ref-lines", "read-lines", "monomer", "monomer-period", "monomer-cut")
FORMATS = ("png", "pdf")


def _get(row, name):
    """Fetch a column, accepting either '-' or '_' in its name; '' means unset."""
    for key in (name, name.replace("-", "_"), name.replace("_", "-")):
        if key in row and row[key] is not None and str(row[key]).strip():
            return str(row[key]).strip()
    return None


def read_manifest(path):
    """Parse and validate. Raises ValueError listing every problem at once."""
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise ValueError(f"{path}: no data rows")
    have = {c.replace("_", "-") for c in rows[0]}
    missing = [c for c in REQUIRED if c.replace("_", "-") not in have]
    if missing:
        raise ValueError(f"{path}: missing required column(s): {', '.join(missing)}")

    out, problems = [], []
    for i, r in enumerate(rows, 1):
        vals = {c: _get(r, c) for c in REQUIRED + OPTIONAL}
        for c in REQUIRED:
            if not vals[c]:
                problems.append(f"row {i}: empty {c}")
        if problems and problems[-1].startswith(f"row {i}:"):
            continue
        for c in ("bam", "reference"):
            if not os.path.exists(vals[c]):
                problems.append(f"row {i} ({vals['suffix']}): {c} not found: {vals[c]}")
        for c in ("k", "min-seg", "merge-gap", "monomer-period"):
            if vals[c] is not None:
                try:
                    int(vals[c])
                except ValueError:
                    problems.append(f"row {i} ({vals['suffix']}): {c} is not an "
                                    f"integer: {vals[c]!r}")
        vals["_row"] = i
        out.append(vals)
    if problems:
        raise ValueError(f"{path}: {len(problems)} unusable row(s):\n  "
                         + "\n  ".join(problems))
    return out


def params_for(v):
    p = Params()
    if v["k"]:
        p.kmer = int(v["k"])
    if v["min-seg"]:
        p.min_seg = int(v["min-seg"])
    if v["merge-gap"]:
        p.merge_gap = int(v["merge-gap"])
    if v["panel-mm"]:
        p.panel_mm = float(v["panel-mm"])
    if v["colour_main"]:
        p.colour_main = v["colour_main"]
    if v["colour_ext"]:
        p.colour_ext = v["colour_ext"]
    if v["monomer"] and v["monomer"].lower() not in ("0", "no", "false", "n"):
        p.monomer = True
    if v["monomer-period"]:
        p.monomer_period = int(v["monomer-period"])
    if v["monomer-cut"]:
        p.monomer_cut = float(v["monomer-cut"])
    return p


def outputs_for(v):
    stem = os.path.join(v["outdir"], v["suffix"])
    return stem, [f"{stem}.quad.{f}" for f in FORMATS]


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="readrefdot-batch", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest", help="TSV, one row per plot")
    ap.add_argument("--force", action="store_true",
                    help="redraw even when the output already exists")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would be drawn, draw nothing")
    a = ap.parse_args(argv)

    try:
        rows = read_manifest(a.manifest)
    except (OSError, ValueError) as e:
        sys.exit(str(e))

    todo, skipped = [], 0
    for v in rows:
        _, paths = outputs_for(v)
        if not a.force and all(os.path.exists(p) for p in paths):
            print(f"  skip   {v['suffix']}  (output exists)")
            skipped += 1
        else:
            todo.append(v)
    if not todo:
        print(f"\nnothing to do: {skipped} row(s) already drawn "
              f"({'--force to redraw' if skipped else ''})")
        return

    # One BAM scan per (bam, reference), not one per row: the scan dominates runtime.
    groups = OrderedDict()
    for v in todo:
        groups.setdefault((v["bam"], v["reference"]), []).append(v)

    n_ok = n_fail = 0
    done = 0
    for (bam, ref), members in groups.items():
        by_read = OrderedDict()
        for v in members:
            by_read.setdefault(v["readid"], []).append(v)
        if a.dry_run:
            for v in members:
                done += 1
                print(f"[{done}/{len(todo)}] would draw {v['suffix']}  "
                      f"({os.path.basename(bam)}, {v['readid']})")
            continue
        for ctx in load(bam, ref, list(by_read), skip_missing=True):
            if isinstance(ctx, ReadNotFound):
                for v in by_read.get(ctx.read_id, []):
                    done += 1
                    print(f"[{done}/{len(todo)}] FAILED {v['suffix']}: {ctx}")
                    n_fail += 1
                continue
            for v in by_read[ctx.read_id]:
                done += 1
                stem, _ = outputs_for(v)
                os.makedirs(v["outdir"], exist_ok=True)
                lines = Lines.parse(v["ref-lines"], v["read-lines"])
                try:
                    paths, st = quad(ctx, params_for(v), stem,
                                     lines=lines or None, formats=FORMATS)
                    print(f"[{done}/{len(todo)}] {v['suffix']}  {ctx.window}  "
                          f"{st['n_fwd']:,} fwd / {st['n_rev']:,} rev  "
                          f"-> {os.path.basename(paths[0])}")
                    n_ok += 1
                except Exception as e:                 # keep going over a long manifest
                    print(f"[{done}/{len(todo)}] FAILED {v['suffix']}: {e}")
                    n_fail += 1

    if not a.dry_run:
        print(f"\n{n_ok} drawn, {skipped} skipped"
              + (f", {n_fail} failed" if n_fail else ""))
        if n_fail:
            sys.exit(1)


if __name__ == "__main__":
    main()
