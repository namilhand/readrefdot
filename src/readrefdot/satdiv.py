"""The divergence plot (--satdiv) - pairwise divergence between the satellite monomers of
a read and its reference window.

A second view of the same data readrefdot already has: the same tiling, the same grouping,
the same annotated intervals, the same Params. `plot.quad` draws the sequence, `draw` here
draws how far apart its monomers are.

The quad dot plot shows where sequence recurs; this shows how far apart the copies are.
Both sequences are cut into CEN178 monomers (the same tiling and the same grouping
readrefdot uses), every monomer is compared with every other, and the matrix of percent
divergence is drawn as a heat map.

The layout is the quad plot's: [reference | read] on BOTH axes, so one matrix fills four
quadrants.

    bottom-left   reference x reference   the reference window's own monomers
    top-right     read x read             the read's own monomers
    top-left      reference (x) x read (y)   every read monomer against every reference one
    bottom-right  read (x) x reference (y)   the transpose

Higher-order repeat structure is what this makes visible. In an array built of a
repeating cassette of m monomers, monomer i and monomer i+m are near-identical while
their neighbours are not, so the matrix carries a ladder of low-divergence lines at
multiples of the HOR period -- off-diagonals the dot plot only hints at.

Divergence is measured in CONSENSUS COLUMNS. Each monomer is aligned to the consensus
and projected onto its 178 columns, so an indel inside one monomer shifts nothing
downstream and any two monomers -- including one from each block -- are compared
base-for-corresponding-base.
"""

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.patches import Rectangle

from . import monomer as mono
from .plot import (AXLAB_PT, COL_AXLAB, MM, STEM_MM, STYLE, _dot_size,
                   _lollipops, block_labels, insertion_sites)

GAP_MM = 0.0               # panel to strip: the sticks start at the edge of the panel
CB_GAP_MM = 3.5            # strip to colour bar
CB_W_MM = 1.6
CHUNK = 64                 # rows of the matrix computed at once
BOX_LW = 0.3               # the box drawn around an annotated interval
COL_ANNOT = "#FFFFFF"      # white: the annotation has to read against a dark heat map
COL_LAB = "#444444"


def project(unit, cons, match=2, mismatch=-3, gap=-5):
    """One monomer written out in the consensus's own columns.

    Returns a uint8 array of len(cons): the unit base aligned to each consensus position,
    or 0 where the unit has lost that position. Bases the unit has inserted are dropped --
    they have no consensus column to sit in, and keeping them would push every later
    position of that one monomer out of register with all the others."""
    got = mono.semiglobal_dp(unit, cons, match, mismatch, gap)
    if got is None:
        return None
    W, C, H = got
    m = len(cons)
    proj = np.zeros(m, dtype=np.uint8)
    i, k = m, int(np.argmax(H[m]))
    while i > 0:
        if k > 0 and H[i, k] == H[i - 1, k - 1] + (match if W[k - 1] == C[i - 1]
                                                   else mismatch):
            proj[i - 1] = W[k - 1]
            i -= 1
            k -= 1
        elif H[i, k] == H[i - 1, k] + gap:
            i -= 1                          # consensus position absent from this monomer
        else:
            k -= 1                          # inserted base: no column for it
    return proj


def divergence(P):
    """Percent divergence between every pair of projected monomers.

    A column counts against the pair when the two differ, including one having lost the
    position and the other not; columns both have lost are not counted either way."""
    n, m = P.shape
    D = np.zeros((n, n), dtype=float)
    if n == 0:
        return D
    isgap = P == 0
    for i0 in range(0, n, CHUNK):
        blk, gb = P[i0:i0 + CHUNK], isgap[i0:i0 + CHUNK]
        both = gb[:, None, :] & isgap[None, :, :]
        diff = ((blk[:, None, :] != P[None, :, :]) & ~both).sum(-1)
        den = m - both.sum(-1)
        D[i0:i0 + CHUNK] = np.where(den > 0, diff / np.maximum(den, 1), 0.0)
    np.fill_diagonal(D, 0.0)
    return D * 100.0


def _colour(u):
    if not u.satellite:
        return mono.NONSAT
    if u.group < 0 or u.group >= len(mono.PALETTE):
        return mono.UNSET
    return mono.PALETTE[u.group]


def _block_units(track, block):
    """The monomers of one block that go into the matrix, in array order.

    Only whole satellite monomers are compared -- the same ones the grouping used. A
    partial unit at the edge of the window is a fragment of a monomer, and a
    non-satellite stretch is not a monomer at all; either would be a spurious row."""
    return [u for u in track.units
            if u.block == block and not u.partial and u.satellite]


def build(S, track, cons):
    """(units, n_ref, divergence matrix) over [reference | read], reference first.

    One matrix for both blocks, so the cross quadrants come out of the same comparison as
    the two self quadrants: a read monomer and a reference monomer are projected onto the
    same consensus columns, so the number in a cross quadrant means exactly what the
    numbers on the diagonal blocks mean."""
    units, P = [], []
    n_ref = 0
    for block in ("ref", "read"):
        for u in _block_units(track, block):
            p = project(S[u.start:u.end], cons)
            if p is None:
                continue
            units.append(u)
            P.append(p)
            n_ref += block == "ref"
    if not units:
        return [], 0, np.zeros((0, 0))
    return units, n_ref, divergence(np.array(P, dtype=np.uint8))


def _index_at(units, pos):
    """Fractional monomer index of a base position: unit 3 half way through is 2.5.

    Interpolating inside the unit rather than snapping to it puts the edge of a drawn box
    on the base the annotation names, not on the nearest monomer boundary."""
    if not units:
        return 0.0
    starts = np.array([u.start for u in units])
    ends = np.array([u.end for u in units])
    if pos <= starts[0]:
        return 0.0
    if pos >= ends[-1]:
        return float(len(units))
    i = int(np.searchsorted(starts, pos, side="right") - 1)
    if pos >= ends[i]:                  # inside a stretch no monomer covers
        return float(i + 1)
    return i + (pos - starts[i]) / (ends[i] - starts[i])


def _spans(ctx, lines, units, n_ref):
    """The annotation in monomer coordinates.

    Returns (reference spans, read spans, insertion marks, deletion marks). Marks are
    monomer coordinates; each belongs to the block it is drawn in.

    `ref-lines`/`read-lines` hold the same intervals readrefdot draws: for an INS the
    donor region on the reference, and in the read both the donor's copy and the inserted
    segment; for a DEL the deleted block on the reference. `Lines.intervals` pairs the
    endpoints (see `annotate.pair_up`) and clips each to its own block.

    An event that is a JUNCTION rather than a segment gets a mark. There are two, and they
    are each other's mirror. A deletion is a junction in the read: it gets no box there,
    and nothing in the divergence matrix says where the sequence was lost -- the monomers
    either side of it are simply neighbours. An insertion is a junction in the reference:
    the read block boxes the inserted copy, but the reference never held it, and the site
    comes from the alignment rather than from the annotation columns."""
    if not lines or not units:
        return [], [], [], []
    ref_iv, read_iv = lines.intervals(ctx)
    ref_u, read_u = units[:n_ref], units[n_ref:]
    out, drawn = [], set()
    for iv, us, off in ((ref_iv, ref_u, 0), (read_iv, read_u, n_ref)):
        got = []
        for a, b in iv:
            if not us or b <= us[0].start or a >= us[-1].end:
                continue                    # outside the monomers of this block
            i, j = _index_at(us, a) + off, _index_at(us, b) + off
            if j - i >= 0.5:                # a point, not a segment (a deletion junction)
                got.append((i, j))
                if off:
                    drawn |= {a, b}
        out.append(got)

    R = len(ctx.ref_seq)
    dels = []
    for p in sorted({R + p for p in lines.read}) if read_u else []:
        if p in drawn or p <= read_u[0].start or p >= read_u[-1].end:
            continue
        dels.append(_index_at(read_u, p) + n_ref)
    ins = []
    for x in insertion_sites(ctx, read_iv) if ref_u else []:
        if ref_u[0].start <= x <= ref_u[-1].end:
            ins.append(_index_at(ref_u, x))
    return out[0], out[1], ins, dels


def _draw_boxes(ax, ref_spans, read_spans):
    """A box around each annotated interval, in every quadrant that holds one.

    An annotated interval -- the donor region, the inserted segment, the deleted block --
    is a run of monomers, and what a run of monomers produces here is a diagonal. The box
    is what that diagonal is FOR. Only the two self quadrants are boxed, each on its own
    diagonal; a cross quadrant would say what those two and the data between them already
    say."""
    rects = [(a, b, a, b) for a, b in ref_spans + read_spans]
    for x0, x1, y0, y1 in rects:
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False,
                               edgecolor=COL_ANNOT, linewidth=BOX_LW, zorder=5))
    return len(rects)


def _draw_marks(ax, ins, dels, n_ref, n):
    """A dashed cross-hair at each junction, inside the block it belongs to: the insertion
    site in the reference, the deletion junction in the read."""
    for marks, lo, hi in ((ins, 0, n_ref), (dels, n_ref, n)):
        for m in marks:
            ax.vlines(m, lo, hi, colors=COL_ANNOT, lw=BOX_LW, ls=(0, (2, 2)), zorder=5)
            ax.hlines(m, lo, hi, colors=COL_ANNOT, lw=BOX_LW, ls=(0, (2, 2)), zorder=5)
    return len(ins) + len(dels)


def scale(params):
    """A stepped sequential scale: one colour per `satdiv_step` %, black above the top.

    Sequential, not diverging: divergence has a true zero and no meaningful midpoint, so a
    diverging map invents a centre and spends half its range on values the data never has.
    viridis is perceptually uniform, so equal steps of divergence look equally different
    (the bands ARE equal steps), and its lightness is monotone, so the ordering survives
    greyscale and any colour vision.

    Fixed rather than fitted to the data, so the same colour means the same divergence in
    every plot and two reads can be compared by eye. A pair further apart than the top of
    the scale takes the top band's own colour: it is off the scale, not off the variable,
    and giving it a colour of its own (black) put the most distant pairs at the dark end,
    where the most alike ones live. Clamping keeps the ramp monotone -- lighter is always
    further apart -- and the title counts how many pairs are held there."""
    vmax, step = params.satdiv_vmax, params.satdiv_step
    bounds = np.arange(0.0, vmax + step / 2, step)
    cmap = plt.get_cmap(params.satdiv_cmap, len(bounds) - 1).copy()
    cmap.set_over(cmap(cmap.N - 1))         # clamp, do not colour differently
    return cmap, BoundaryNorm(bounds, cmap.N)


def _strips(fig, geom, units, n_ref, n, panel_mm):
    """The monomer annotation: one lollipop per monomer, along the top and the right of
    the whole panel, so it labels the reference block and the read block in turn."""
    x0, y0, box, gap, strip, fw, fh = geom
    centres = np.arange(n) + 0.5
    cols = [_colour(u) for u in units]
    dot = _dot_size(panel_mm, n)
    stem = STEM_MM * MM / strip
    top = fig.add_axes([x0 / fw, (y0 + box + gap) / fh, box / fw, strip / fh])
    right = fig.add_axes([(x0 + box + gap) / fw, y0 / fh, strip / fw, box / fh])
    _lollipops(top, centres, cols, dot, True, stem)
    _lollipops(right, centres, cols, dot, False, stem)
    top.set_xlim(0, n); top.set_ylim(0, 1)
    right.set_xlim(0, 1); right.set_ylim(0, n)
    top.axvline(n_ref, color="black", lw=0.5, ymax=stem)
    right.axhline(n_ref, color="black", lw=0.5, xmax=stem)
    for a in (top, right):
        a.set_xticks([]); a.set_yticks([])
        for sp in a.spines.values():        # the panel's own frame is the strip's base
            sp.set_visible(False)
    return top


def draw(ctx, params, out_stem, lines=None, formats=("pdf", "png")):
    """Draw the quad divergence plot. Returns (paths, stats)."""
    S = ctx.ref_seq + ctx.read_seq
    track = mono.annotate(ctx, period=params.monomer_period, cut=params.monomer_cut,
                          consensus=params.monomer_consensus)
    if track is None:
        raise ValueError("not a tandem satellite array: no monomers to compare")
    cons = params.monomer_consensus
    if cons:
        # annotate() may have phased the array on the reverse strand; project onto the
        # same orientation it tiled with, or every unit would align back to front.
        cons = mono.revcomp(cons) if "reverse" in track.phase else cons
    else:
        cons = S[track.units[0].start:track.units[0].end]

    units, n_ref, D = build(S, track, cons)
    n = len(units)
    if n == 0:
        raise ValueError("no whole satellite monomers to compare")
    cmap, norm = scale(params)

    mpl.rcParams.update(STYLE)
    panel_mm = params.satdiv_panel_mm or params.panel_mm
    box = panel_mm * MM
    dot = _dot_size(panel_mm, n)
    strip = STEM_MM * MM + dot / 72 / 2 + 0.05 * MM     # the stick plus half a circle
    gap, cb_gap, cb_w = GAP_MM * MM, CB_GAP_MM * MM, CB_W_MM * MM
    ml, mb, mt, mr = 0.17, 0.17, 0.50, 0.34
    fw = ml + box + gap + strip + cb_gap + cb_w + mr
    fh = mb + box + gap + strip + mt
    fig = plt.figure(figsize=(fw, fh))

    ax = fig.add_axes([ml / fw, mb / fh, box / fw, box / fh])
    # interpolation="none" embeds the matrix at its own n x n size in the vector output
    # instead of a resampled copy at device resolution: one PDF sample per monomer pair,
    # which is both smaller and exactly the data.
    im = ax.imshow(D, cmap=cmap, norm=norm, extent=(0, n, 0, n), origin="lower",
                   interpolation="none", aspect="auto")
    ref_spans, read_spans, ins, dels = _spans(ctx, lines, units, n_ref)
    n_box = _draw_boxes(ax, ref_spans, read_spans)
    n_mark = _draw_marks(ax, ins, dels, n_ref, n)
    ax.set_xlim(0, n); ax.set_ylim(0, n)
    ax.axvline(n_ref, color=COL_ANNOT, lw=0.5, zorder=4)
    ax.axhline(n_ref, color=COL_ANNOT, lw=0.5, zorder=4)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.5); sp.set_color("black")

    strip_ax = _strips(fig, (ml, mb, box, gap, strip, fw, fh), units, n_ref, n,
                       panel_mm)

    for pos, lab in zip((n_ref / 2, n_ref + (n - n_ref) / 2), block_labels(ctx)):
        ax.annotate(lab, xy=(pos, 0), xycoords=("data", "axes fraction"),
                    xytext=(0, -6), textcoords="offset points",
                    ha="center", va="top", fontsize=AXLAB_PT, color=COL_AXLAB)
        ax.annotate(lab, xy=(0, pos), xycoords=("axes fraction", "data"),
                    xytext=(-8, 0), textcoords="offset points", rotation=90,
                    ha="right", va="center", fontsize=AXLAB_PT, color=COL_AXLAB)

    cax = fig.add_axes([(ml + box + gap + strip + cb_gap) / fw, mb / fh, cb_w / fw,
                        box / fh])
    tick = next(t for t in (1, 2, 4, 5, 10, 20, 50) if params.satdiv_vmax / t <= 6)
    ticks = np.arange(0, params.satdiv_vmax + tick / 2, tick)
    cb = fig.colorbar(im, cax=cax, ticks=ticks)
    # The top band holds everything above it, so say so on the bar rather than with a
    # colour of its own.
    cb.set_ticklabels([f"{t:g}" for t in ticks[:-1]] + [f"\u2265{ticks[-1]:g}"])
    cb.outline.set_linewidth(0.4)
    cax.tick_params(width=0.4, length=1.6, labelsize=4.5, pad=1.2)
    cax.set_ylabel("monomer-pair divergence (%)", fontsize=5, color=COL_LAB, labelpad=2)

    over = int((D > params.satdiv_vmax).sum())
    title = (f"{ctx.read_id}\n{ctx.window}  (strand {ctx.strand})\n"
             f"{n_ref} + {n - n_ref} monomers of {track.period} bp, "
             f"{track.n_groups} groups at {int(params.monomer_cut * 100)}% identity"
             f"\nphase: {track.phase}"
             + (f"  ·  {track.n_nonsatellite} non-satellite" if track.n_nonsatellite
                else "")
             + (f"  ·  {over:,} pairs over {params.satdiv_vmax:g}%" if over else ""))
    if n_box or n_mark:
        title += ("\nbox: annotated donor, inserted or deleted segment"
                  if n_box else "")
        notes = (["the insertion site in the reference"] if ins else []) + \
                (["the deletion junction in the read"] if dels else [])
        title += ("\ndashed: " + " and ".join(notes)) if notes else ""
    strip_ax.set_title(title, fontsize=5, linespacing=1.6)

    paths = []
    for fmt in formats:
        path = f"{out_stem}.satdiv.{fmt}"
        # With pdf.compression on, matplotlib turns any image of 256 colours or fewer
        # into a 4-bit INDEXED-palette image. A stepped scale has about twenty colours,
        # so this plot always trips it -- and Illustrator drops the palette when the PDF
        # is placed, which is why the heat map arrived there with no colour at all. The
        # dot plot is pure vector and never had the problem. Writing this PDF
        # uncompressed keeps the image in plain DeviceRGB; it costs ~0.2 MB.
        with mpl.rc_context({"pdf.compression": 0} if fmt == "pdf" else {}):
            fig.savefig(path, format=fmt, dpi=params.dpi)
        paths.append(path)
    plt.close(fig)
    return paths, dict(n_ref=n_ref, n_read=n - n_ref, n_boxes=n_box, n_marks=n_mark,
                       n_over=over,
                       monomer=track, units=units, matrix=D)


def write_tsv(ctx, units, n_ref, D, path):
    """The matrix itself: one row per monomer over [reference | read], in plot order."""
    R = len(ctx.ref_seq)
    with open(path, "w") as fh:
        head = "\t".join(str(i + 1) for i in range(len(units)))
        fh.write(f"block\tindex\tblock_start\tblock_end\tgroup\t{head}\n")
        for i, u in enumerate(units):
            b = u.start if i < n_ref else u.start - R
            vals = "\t".join(f"{v:.3f}" for v in D[i])
            fh.write(f"{u.block}\t{i + 1}\t{b}\t{b + u.length}\t"
                     f"{u.group if u.group >= 0 else '.'}\t{vals}\n")
    return path
