"""readrefdot - k-mer dot plots of a long read against the reference it maps to."""

import argparse
import os
import re
import sys

from . import __version__
from .annotate import Lines
from . import tree as tree_mod
from .monomer import CEN178, DEFAULT_CUT, write_tsv
from .plot import Params, quad
from .read import ReadNotFound, iter_primary, load

_SAFE = re.compile(r"[^A-Za-z0-9._-]+")


def safe_name(read_id):
    """Read ids contain '/', which cannot go in a filename."""
    return _SAFE.sub("_", read_id).strip("_")


def read_consensus(arg):
    """None -> the built-in CEN178; 'none' -> de novo phase; else the first FASTA record."""
    if arg is None:
        return CEN178
    if arg.lower() in ("none", "off", "-"):
        return None
    seq = []
    with open(arg) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq:
                    break
                continue
            seq.append(line.strip())
    text = "".join(seq).upper()
    if not text:
        sys.exit(f"--monomer-consensus: no sequence in {arg}")
    return text


def build_parser():
    p = argparse.ArgumentParser(
        prog="readrefdot",
        description="k-mer dot plots comparing a long read with the reference genome "
                    "region it maps to.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--version", action="version", version=f"readrefdot {__version__}")
    p.add_argument("--bam", required=True, help="BAM containing the read(s)")
    p.add_argument("--ref", required=True, help="reference FASTA the BAM was aligned to")

    sel = p.add_mutually_exclusive_group(required=True)
    sel.add_argument("--read", nargs="+", metavar="ID", help="read id(s) to plot")
    sel.add_argument("--all", action="store_true",
                     help="plot every primary alignment in the BAM")

    p.add_argument("--outdir", default=".", help="directory for the output files")
    p.add_argument("--name", metavar="NAME",
                   help="output file stem, giving <NAME>.quad.png/.pdf "
                        "(default: the read id)")
    p.add_argument("-k", "--kmer", type=int, default=20, help="k-mer size (5-31)")
    p.add_argument("--min-seg", type=int, default=170,
                   help="drop diagonal runs shorter than this (bp)")
    p.add_argument("--merge-gap", type=int, default=None,
                   help="max gap chained into one run (default k+1, which bridges a "
                        "single substitution)")
    p.add_argument("--panel-mm", "--panel_mm", type=float, default=45.0,
                   help="size of the square plot box in mm, excluding title and labels")
    p.add_argument("--colour_main", "--colour-main", default=None, metavar="COLOUR",
                   help="colour of the diagonals the aligner placed the read on "
                        "(default: black)")
    p.add_argument("--colour_ext", "--colour-ext", default=None, metavar="COLOUR",
                   help="colour of every other diagonal (default: grey40)")
    p.add_argument("--monomer", action="store_true",
                   help="annotate satellite (CEN178) monomers along the top and right "
                        "of the panel instead of genomic coordinates")
    p.add_argument("--monomer-period", type=int, default=None, metavar="BP",
                   help="satellite unit length (default: detect it from the sequence)")
    p.add_argument("--monomer-style", choices=("lollipop", "block"), default="lollipop",
                   help="how a monomer is drawn on the axis")
    p.add_argument("--monomer-consensus", default=None, metavar="FASTA",
                   help="repeat consensus that fixes where a unit starts "
                        "(default: the published CEN178 monomer; 'none' = take the phase "
                        "from the sequence itself)")
    p.add_argument("--tree-method", choices=("dendrogram", "nj"), default="dendrogram",
                   help="dendrogram = the grouping tree itself (colours are its branches); "
                        "nj = a separate neighbour-joining tree")
    p.add_argument("--no-tree", "--no-dendrogram", action="store_true", dest="no_tree",
                   help="skip the monomer dendrogram")
    p.add_argument("--monomer-cut", type=float, default=DEFAULT_CUT, metavar="F",
                   help="group monomers whose estimated identity is at least this")
    p.add_argument("--monomer-tsv", action="store_true",
                   help="also write <NAME>.monomers.tsv, one row per unit")
    p.add_argument("--annot-style", choices=("box", "lines", "both"), default="box",
                   help="how --ref-lines/--read-lines are drawn: a black box around each "
                        "annotated interval, dotted guide lines at its ends, or both")
    p.add_argument("--dpi", type=int, default=None,
                   help="resolution of the PNG; the PDF stays vector either way "
                        "(default: 600)")
    p.add_argument("--ref-lines", metavar="P,...",
                   help="annotate the interval(s) between these reference "
                        "positions (1-based)")
    p.add_argument("--read-lines", metavar="P,...",
                   help="annotate the interval(s) between these read "
                        "positions (0-based)")
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    if not 5 <= a.kmer <= 31:
        sys.exit("--kmer must be between 5 and 31")
    lines = Lines.parse(a.ref_lines, a.read_lines)
    if lines and a.all:
        sys.exit("--ref-lines/--read-lines describe one read; use --read, not --all")
    if lines and len(a.read) > 1:
        sys.exit("--ref-lines/--read-lines describe one read; give a single --read")
    if a.name and (a.all or len(a.read) > 1):
        sys.exit("--name gives one file stem; use it with a single --read")

    read_ids = list(iter_primary(a.bam)) if a.all else a.read
    if not read_ids:
        sys.exit(f"no primary alignments in {a.bam}")
    os.makedirs(a.outdir, exist_ok=True)
    params = Params(kmer=a.kmer, min_seg=a.min_seg, merge_gap=a.merge_gap,
                    panel_mm=a.panel_mm, monomer=a.monomer or a.monomer_tsv,
                    monomer_period=a.monomer_period, monomer_cut=a.monomer_cut,
                    monomer_style=a.monomer_style,
                    monomer_consensus=read_consensus(a.monomer_consensus),
                    annot_style=a.annot_style)
    if a.dpi:
        params.dpi = a.dpi
    if a.colour_main:
        params.colour_main = a.colour_main
    if a.colour_ext:
        params.colour_ext = a.colour_ext

    n_ok = n_fail = 0
    for i, ctx in enumerate(load(a.bam, a.ref, read_ids, skip_missing=True), 1):
        if isinstance(ctx, ReadNotFound):
            print(f"[{i}/{len(read_ids)}] FAILED: {ctx}", file=sys.stderr)
            n_fail += 1
            continue
        stem = os.path.join(a.outdir, a.name or safe_name(ctx.read_id))
        paths, st = quad(ctx, params, stem, lines=lines or None)
        if st["monomer"] is not None and not a.no_tree:
            tree_mod.draw(st["monomer"], stem, panel_mm=a.panel_mm,
                          method=a.tree_method,
                          title=f"{ctx.read_id}\n{st['monomer'].n_full} monomers, "
                                f"{st['monomer'].n_groups} groups")
        if a.monomer_tsv and st["monomer"] is not None:
            write_tsv(st["monomer"], ctx, f"{stem}.monomers.tsv")
        if params.monomer and st["monomer"] is None:
            print(f"    note: {ctx.read_id} is not a tandem array; "
                  f"drew coordinates instead", file=sys.stderr)
        print(f"[{i}/{len(read_ids)}] {ctx.read_id}  {ctx.window}  "
              f"ref {len(ctx.ref_seq):,} + read {ctx.read_len:,} bp  "
              f"{st['n_fwd']:,} fwd / {st['n_rev']:,} rev  -> {os.path.basename(paths[0])}")
        n_ok += 1

    print(f"\n{n_ok} plot(s) in {a.outdir}" + (f", {n_fail} failed" if n_fail else ""))
    if n_fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
