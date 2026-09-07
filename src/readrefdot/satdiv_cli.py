"""satdivplot - pairwise monomer divergence heat maps, one page per read.

    satdivplot --bam s.bam --ref genome.fa --read ID [--outdir DIR]
    satdivplot-batch manifest.tsv

The manifest is the same file readrefdot-batch takes; the columns it uses are
bam, readid, reference, outdir, suffix, monomer-period, monomer-cut, monomer-consensus,
plus satdiv-panel-mm and satdiv-cmap for this plot's own settings.
"""

import argparse
import os
import sys
from collections import OrderedDict

from . import __version__
from . import satdiv
from .batch import read_manifest
from .cli import read_consensus, safe_name
from .read import ReadNotFound, iter_primary, load

FORMATS = ("pdf", "png")


def _add_common(p):
    p.add_argument("--panel-mm", "--panel_mm", type=float, default=satdiv.PANEL_MM,
                   help="size of each square heat map in mm")
    p.add_argument("--cmap", default=satdiv.CMAP,
                   help="matplotlib colormap for divergence (diverging by default)")
    p.add_argument("--vmin", type=float, default=None,
                   help="%% divergence at the low end of the colour scale "
                        "(default: the 2nd percentile of the data)")
    p.add_argument("--vmax", type=float, default=None,
                   help="%% divergence at the high end (default: the 98th percentile)")
    p.add_argument("--center", "--centre", type=float, default=None, dest="center",
                   help="%% divergence the diverging scale is neutral at "
                        "(default: the array's median)")
    p.add_argument("--monomer-period", type=int, default=None, metavar="BP",
                   help="satellite unit length (default: detect it from the sequence)")
    p.add_argument("--monomer-cut", type=float, default=satdiv.mono.DEFAULT_CUT,
                   metavar="F", help="group monomers whose estimated identity is at "
                                     "least this")
    p.add_argument("--monomer-consensus", default=None, metavar="FASTA",
                   help="repeat consensus that fixes where a unit starts, and the "
                        "columns divergence is measured in "
                        "(default: the published CEN178 monomer)")
    p.add_argument("--matrix-tsv", action="store_true",
                   help="also write <NAME>.satdiv.{ref,read}.tsv, the matrices themselves")


def _params(a):
    return satdiv.Params(panel_mm=a.panel_mm, cmap=a.cmap, vmin=a.vmin, vmax=a.vmax,
                         center=a.center,
                         monomer_period=a.monomer_period, monomer_cut=a.monomer_cut,
                         monomer_consensus=read_consensus(a.monomer_consensus))


def _draw(ctx, params, stem, matrix_tsv=False):
    paths, st = satdiv.draw(ctx, params, stem, formats=FORMATS)
    if matrix_tsv:
        for name in ("ref", "read"):
            units = [u for u in st["monomer"].units
                     if u.block == name and not u.partial and u.satellite]
            satdiv.write_tsv(ctx, st["monomer"], name, units, st["matrices"][name],
                             f"{stem}.satdiv.{name}.tsv")
    return paths, st


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="satdivplot", description=__doc__.split("\n\n")[0],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--version", action="version", version=f"satdivplot {__version__}")
    p.add_argument("--bam", required=True, help="BAM containing the read(s)")
    p.add_argument("--ref", required=True, help="reference FASTA the BAM was aligned to")
    sel = p.add_mutually_exclusive_group(required=True)
    sel.add_argument("--read", nargs="+", metavar="ID", help="read id(s) to plot")
    sel.add_argument("--all", action="store_true",
                     help="plot every primary alignment in the BAM")
    p.add_argument("--outdir", default=".", help="directory for the output files")
    p.add_argument("--name", metavar="NAME",
                   help="output file stem (default: the read id)")
    _add_common(p)
    a = p.parse_args(argv)
    if a.name and (a.all or len(a.read) > 1):
        sys.exit("--name gives one file stem; use it with a single --read")

    read_ids = list(iter_primary(a.bam)) if a.all else a.read
    if not read_ids:
        sys.exit(f"no primary alignments in {a.bam}")
    os.makedirs(a.outdir, exist_ok=True)
    params = _params(a)

    n_ok = n_fail = 0
    for i, ctx in enumerate(load(a.bam, a.ref, read_ids, skip_missing=True), 1):
        if isinstance(ctx, ReadNotFound):
            print(f"[{i}/{len(read_ids)}] FAILED: {ctx}", file=sys.stderr)
            n_fail += 1
            continue
        stem = os.path.join(a.outdir, a.name or safe_name(ctx.read_id))
        try:
            paths, st = _draw(ctx, params, stem, a.matrix_tsv)
        except Exception as e:
            print(f"[{i}/{len(read_ids)}] FAILED {ctx.read_id}: {e}", file=sys.stderr)
            n_fail += 1
            continue
        print(f"[{i}/{len(read_ids)}] {ctx.read_id}  {ctx.window}  "
              f"ref {st['n_ref']} + read {st['n_read']} monomers  "
              f"{st['vmin']:.1f}-{st['vmax']:.1f}% scale  "
              f"-> {os.path.basename(paths[0])}")
        n_ok += 1

    print(f"\n{n_ok} plot(s) in {a.outdir}" + (f", {n_fail} failed" if n_fail else ""))
    if n_fail:
        sys.exit(1)


def _params_for(v):
    p = satdiv.Params()
    if v["satdiv-panel-mm"]:
        p.panel_mm = float(v["satdiv-panel-mm"])
    if v["satdiv-cmap"]:
        p.cmap = v["satdiv-cmap"]
    if v["monomer-period"]:
        p.monomer_period = int(v["monomer-period"])
    if v["monomer-cut"]:
        p.monomer_cut = float(v["monomer-cut"])
    if v["monomer-consensus"]:
        p.monomer_consensus = read_consensus(v["monomer-consensus"])
    return p


def _outputs(v):
    stem = os.path.join(v["outdir"], v["suffix"])
    return stem, [f"{stem}.satdiv.{f}" for f in FORMATS]


def batch_main(argv=None):
    ap = argparse.ArgumentParser(
        prog="satdivplot-batch", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest", help="TSV, one row per plot")
    ap.add_argument("--outdir", default=None,
                    help="write everything here instead of the manifest's outdir column")
    ap.add_argument("--matrix-tsv", action="store_true",
                    help="also write the matrices as TSV beside each plot")
    ap.add_argument("--force", action="store_true",
                    help="redraw even when the output already exists")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would be drawn, draw nothing")
    a = ap.parse_args(argv)

    try:
        rows = read_manifest(a.manifest)
    except (OSError, ValueError) as e:
        sys.exit(str(e))
    if a.outdir:
        for v in rows:
            v["outdir"] = a.outdir

    todo, skipped = [], 0
    for v in rows:
        _, paths = _outputs(v)
        if not a.force and all(os.path.exists(p) for p in paths):
            print(f"  skip   {v['suffix']}  (output exists)")
            skipped += 1
        else:
            todo.append(v)
    if not todo:
        print(f"\nnothing to do: {skipped} row(s) already drawn "
              f"({'--force to redraw' if skipped else ''})")
        return

    groups = OrderedDict()
    for v in todo:
        groups.setdefault((v["bam"], v["reference"]), []).append(v)

    n_ok = n_fail = done = 0
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
                stem, _ = _outputs(v)
                os.makedirs(v["outdir"], exist_ok=True)
                try:
                    paths, st = _draw(ctx, _params_for(v), stem, a.matrix_tsv)
                    print(f"[{done}/{len(todo)}] {v['suffix']}  {ctx.window}  "
                          f"ref {st['n_ref']} + read {st['n_read']} monomers  "
                          f"-> {os.path.basename(paths[0])}")
                    n_ok += 1
                except Exception as e:              # keep going over a long manifest
                    print(f"[{done}/{len(todo)}] FAILED {v['suffix']}: {e}")
                    n_fail += 1

    if not a.dry_run:
        print(f"\n{n_ok} drawn, {skipped} skipped" + (f", {n_fail} failed" if n_fail
                                                      else ""))
        if n_fail:
            sys.exit(1)


if __name__ == "__main__":
    main()
