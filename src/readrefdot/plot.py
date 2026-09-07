"""The quad dot plot.

[reference | read] on BOTH axes, from a single self-comparison of the two concatenated
sequences. That one comparison fills four quadrants:

    bottom-left   reference x reference   the reference window's own repeat structure
    top-right     read x read             the read's own repeat structure
    top-left      reference (x) x read (y)   read vs reference
    bottom-right  read (x) x reference (y)   the transpose

Everything drawn is an exact k-mer match found in the sequences themselves. The
aligner's decisions are not shown.
"""

from dataclasses import dataclass

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Rectangle
from matplotlib.transforms import ScaledTranslation

from . import monomer as mono
from .kmer import (build_index, filter_min_length, kmer_codes, match,
                   merge_segments, revcomp)

MM = 1 / 25.4
PANEL_MM = 45.0                 # default plot box, excluding title and labels
TICK_STEP = 5000
STRIP_BLOCK_MM = 1.5            # thickness of the monomer strip, block style
STRIP_GAP_MM = 0.0              # the sticks start at the edge of the panel
DOT_MIN, DOT_MAX = 0.7, 3.0     # pt: monomer marker diameter, auto-sized to the spacing
STEM_MM = 0.66                  # length of a lollipop stick; the circle sits on its end

COL_MAIN = "#000000"            # diagonals the aligner placed the read on
COL_EXT = "#B2B2B2"             # grey70 - every other diagonal
COL_REV = "#D55E00"             # vermillion - reverse-complement matches
COL_LINE = "#0072B2"            # blue - annotation guide lines
COL_BOX = "#0073b2"             # blue - the box around an annotated interval
BOX_LW = 0.3
DPI = 600                       # 45 mm of dot plot is finer than 300 dpi resolves
COL_LAB = "#444444"
COL_AXLAB = "#000000"           # the two block labels under and beside the panel
AXLAB_PT = 5
BASE_PT = 5                     # the text size the millimetre geometry is drawn for
PNG_PANEL_MM = 140.0            # the PNG is the same figure at poster size,
PNG_TEXT_PT = 18.0              # with text this big
DIAG_TOL = 3                    # bp slack when matching a run to an alignment diagonal
MIN_BLOCK = 20                  # smallest aligned block that contributes an identity band
MIN_OVERLAP = 0.5               # a run must lie this far inside a band to count as main
EXTEND_GAP = 200                # bp: bands grow along their diagonal through runs this close

STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7,
    "xtick.labelsize": 5, "ytick.labelsize": 5,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 3, "ytick.major.size": 3,
    "xtick.direction": "out", "ytick.direction": "out",
    "axes.grid": False,
    "savefig.dpi": 300, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42, "ps.fonttype": 42,
}


class Sizes:
    """Every size in one drawing of a figure.

    The PDF is drawn at `--panel-mm` with BASE_PT text, the PNG at PNG_PANEL_MM with
    PNG_TEXT_PT -- which is not the same zoom. Text is deliberately larger relative to the
    panel there, because a PNG is looked at on a screen or a slide rather than placed at
    figure size in a page layout.

    `g` scales anything geometric: a length in mm or inches, a line width, a marker
    diameter. `t` scales anything textual: a font size, a text offset, a pad."""

    def __init__(self, panel_mm, text_pt=BASE_PT):
        self.panel_mm = float(panel_mm)
        self.g = self.panel_mm / PANEL_MM
        self.t = float(text_pt) / BASE_PT

    def rc(self):
        """STYLE with its sizes scaled, for `mpl.rcParams.update`."""
        rc = dict(STYLE)
        for key in ("font.size", "axes.labelsize", "xtick.labelsize", "ytick.labelsize"):
            rc[key] = STYLE[key] * self.t
        for key in ("axes.linewidth", "xtick.major.width", "ytick.major.width",
                    "xtick.major.size", "ytick.major.size"):
            rc[key] = STYLE[key] * self.g
        return rc


def sizes_for(fmt, params, panel_mm=None):
    """The sizes one output format is drawn at."""
    if fmt == "png":
        return Sizes(params.png_panel_mm, params.png_text_pt)
    return Sizes(panel_mm if panel_mm else params.panel_mm)


@dataclass
class Params:
    kmer: int = 20
    min_seg: int = 170          # ~one CEN178 satellite unit (units vary 177-179 bp)
    merge_gap: int = None       # default k+1: bridges a single substitution
    colour_main: str = COL_MAIN
    colour_ext: str = COL_EXT
    panel_mm: float = PANEL_MM
    monomer: bool = False       # annotate satellite monomers instead of coordinates
    monomer_period: int = None  # unit length in bp (default: detect it)
    monomer_cut: float = mono.DEFAULT_CUT   # identity at which units group together
    monomer_style: str = "lollipop"         # "lollipop" or "block"
    monomer_consensus: str = mono.CEN178    # phase reference; None = take it from the data
    annot_style: str = "box"    # "box", "lines" or "both": how --ref/read-lines are drawn
    dpi: int = DPI              # raster resolution; the PDF stays vector either way
    satdiv: bool = False        # also draw the monomer divergence plot (satdiv.draw)
    satdiv_panel_mm: float = None   # its plot box; None = the same as this one
    satdiv_cmap: str = "viridis"    # sequential: dark = alike, light = far apart
    satdiv_vmax: float = 16.0       # % divergence at the top of the scale
    satdiv_step: float = 1.0        # % divergence per colour band
    png_panel_mm: float = PNG_PANEL_MM      # the PNG copy is drawn at this size
    png_text_pt: float = PNG_TEXT_PT        # with text this big

    @property
    def gap(self):
        return self.merge_gap if self.merge_gap is not None else self.kmer + 1


def _seq_x(q, aln, ctx):
    """CIGAR query coord (that alignment's orientation, hard clips included) -> index
    into the plotted read_seq."""
    if aln.is_reverse != ctx.primary.is_reverse:
        q = ctx.full_read_len - q
    return q - ctx.read_offset


def identity_bands(ctx, R, min_block=MIN_BLOCK):
    """Where the ALIGNER placed the read, as (diagonal, x_lo, x_hi) bands.

    A run counts as 'main' only if it lies on one of these diagonals AND overlaps that
    band's x-range. Testing the diagonal alone is not enough: in a repeat array a ladder
    run can sit at the same diagonal offset while being nowhere near the alignment.

    Bands are then grown along their own diagonal through runs that are collinear and
    nearly contiguous (see extend_bands), so a black line is a literally continuous
    diagonal rather than only the aligned part of one. That is the point: when an
    insertion duplicates the reference immediately upstream, the duplicated copy lies on
    the same diagonal, and drawing it black shows how much of the reference the read
    carries twice.

    This is the only use made of the aligner's answer -- it colours runs, it is not
    drawn."""
    n = R + len(ctx.read_seq)
    bands = [(0, 0, n)]                       # self-identity, spans the whole plot
    for a in ctx.alignments:
        if a.is_reverse != ctx.primary.is_reverse:
            continue                          # opposite orientation -> reverse matches
        q, r = 0, a.reference_start
        for op, ln in (a.cigartuples or []):
            if op in (0, 7, 8):
                if ln >= min_block:
                    qs = _seq_x(q, a, ctx)
                    rs = r - ctx.win_start
                    d = (R + qs) - rs
                    bands.append((d, rs, rs + ln))            # x on the reference side
                    bands.append((-d, R + qs, R + qs + ln))   # x on the read side
                q += ln; r += ln
            elif op == 1:
                q += ln
            elif op in (2, 3):
                r += ln
            elif op in (4, 5):
                q += ln
    return bands


def _self_compare(ctx, p):
    """Exact k-mer runs within [reference + read], junction-spanning k-mers excluded."""
    k = p.kmer
    R = len(ctx.ref_seq)
    S = ctx.ref_seq + ctx.read_seq
    n = len(S)

    vals, valid = kmer_codes(S, k)
    valid[max(0, R - k + 1):R] = False
    _, _, sorted_vals, order = build_index(S, k, valid)

    fa, fb, _ = match(vals, valid, sorted_vals, order)
    cv, cvalid = kmer_codes(revcomp(S), k)
    bad = n - k - np.arange(max(0, R - k + 1), R)
    cvalid[bad[(bad >= 0) & (bad < cvalid.size)]] = False
    cq, rb, _ = match(cv, cvalid, sorted_vals, order)
    ra = (n - k - cq).astype(np.int64)

    fwd = filter_min_length(merge_segments(fa, fb, p.gap, False), k, p.min_seg)
    rev = filter_min_length(merge_segments(ra, rb, p.gap, True), k, p.min_seg)
    return fwd, rev, n


def extend_bands(bands, a0, x_end, d, tol=DIAG_TOL, gap=EXTEND_GAP):
    """Grow each identity band along its own diagonal through collinear runs separated
    by at most `gap`, so the black path follows the whole continuous diagonal."""
    out = []
    for bd, lo, hi in bands:
        on = np.abs(d - bd) <= tol
        if on.any():
            xs, xe = a0[on], x_end[on]
            while True:
                near = (xe >= lo - gap) & (xs <= hi + gap)
                if not near.any():
                    break
                nlo, nhi = min(lo, int(xs[near].min())), max(hi, int(xe[near].max()))
                if nlo == lo and nhi == hi:
                    break
                lo, hi = nlo, nhi
        out.append((bd, lo, hi))
    return out


def _block_ticks(plot_lo, plot_hi, coord_lo, fmt):
    """Ticks every TICK_STEP anchored on the coordinate system, labelled only at the
    two ends -- at 45 mm there is no room for more, and the ends give the range."""
    span = plot_hi - plot_lo
    first = -(-coord_lo // TICK_STEP) * TICK_STEP          # round up to a 5 kb multiple
    inner = [plot_lo + (c - coord_lo) for c in
             range(int(first), int(coord_lo + span) + 1, TICK_STEP)]
    pos = [plot_lo] + [p for p in inner if plot_lo < p < plot_hi] + [plot_hi]
    lab = [fmt(coord_lo)] + [""] * (len(pos) - 2) + [fmt(coord_lo + span)]
    return pos, lab


def _mb(v):
    """x.xxx Mb, rounded to the nearest kb (6,981,972 -> 6.982)."""
    return f"{v/1e6:.3f}"


def _ref_ticks(ctx, R):
    """Ticks every 5 kb anchored on the START OF THE ALIGNMENT, not on the window edge,
    so the first tick marks where the read's alignment begins."""
    pos, lab = [], []
    c = ctx.aln_start
    while c <= ctx.win_end:
        p = c - ctx.win_start
        if 0 <= p <= R:
            pos.append(p); lab.append("")
        c += TICK_STEP
    if pos:
        lab[0] = _mb(ctx.aln_start)
        lab[-1] = _mb(ctx.aln_start + (len(pos) - 1) * TICK_STEP)
    return pos, lab


def _ticks(ctx, R, Q):
    """Reference block in Mb (x.xxx, floored), read block in kb, both stepped by 5 kb.
    Returns (positions, labels, n_ref) -- n_ref splits the two blocks so the labelled
    ends can be aligned away from each other."""
    rp, rl = _ref_ticks(ctx, R)
    qp, ql = _block_ticks(R, R + Q, ctx.read_offset, lambda v: f"{v/1000:.1f}")
    return rp + qp, rl + ql, len(rp)


def _dot_size(S, n_units):
    """Marker diameter in pt, from the space one unit actually gets. A 45 mm panel over
    ~220 monomers gives each 0.2 mm, so the dots are sized to sit side by side rather
    than to a fixed size that would overlap into a solid bar."""
    return float(np.clip(S.panel_mm / max(n_units, 1) * 0.9 * 72 / 25.4,
                         DOT_MIN * S.g, DOT_MAX * S.g))


def _lollipops(ax, centres, cols, dot, vertical, tip, lw=0.2):
    """A stick from the panel edge out to a coloured circle centred on the stick's end.

    `tip` is the stick's length as a fraction of the strip."""
    if vertical:                                    # the strip along the top
        segs = [[(c, 0.0), (c, tip)] for c in centres]
        pts = (centres, np.full(centres.size, tip))
    else:                                           # the strip along the right
        segs = [[(0.0, c), (tip, c)] for c in centres]
        pts = (np.full(centres.size, tip), centres)
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=lw, zorder=1))
    ax.scatter(pts[0], pts[1], s=dot ** 2, c=cols, linewidths=0, zorder=2)


def _strips(fig, geom, R, n, track, style, S):
    """The monomer annotation: one mark per satellite unit, along the top and the right
    of the whole panel (so it labels the reference block and the read block in turn).

    Colour = similarity group, so a repeating colour pattern is the array's higher-order
    structure, and the two blocks can be read against each other."""
    ml, mb, box, gap, strip, fw, fh = geom
    starts = np.array([u.start for u in track.units], dtype=float)
    widths = np.array([u.length for u in track.units], dtype=float)
    centres = starts + widths / 2
    cols = track.colours()
    dot = _dot_size(S, len(track.units))
    stem = STEM_MM * MM * S.g / strip               # stick length, as a strip fraction

    top = fig.add_axes([ml / fw, (mb + box + gap) / fh, box / fw, strip / fh])
    right = fig.add_axes([(ml + box + gap) / fw, mb / fh, strip / fw, box / fh])
    if style == "block":
        top.bar(centres, height=1.0, width=widths, color=cols, linewidth=0,
                align="center")
        right.barh(centres, width=1.0, height=widths, color=cols, linewidth=0,
                   align="center")
    else:
        _lollipops(top, centres, cols, dot, True, stem, 0.2 * S.g)
        _lollipops(right, centres, cols, dot, False, stem, 0.2 * S.g)
    top.set_xlim(0, n); top.set_ylim(0, 1)
    right.set_xlim(0, 1); right.set_ylim(0, n)
    top.axvline(R, color="black", lw=0.5 * S.g, ymax=1.0 if style == "block" else stem)
    right.axhline(R, color="black", lw=0.5 * S.g, xmax=1.0 if style == "block" else stem)

    for a in (top, right):
        a.set_xticks([]); a.set_yticks([])
        for sp in a.spines.values():        # the panel's own frame is the strip's base
            sp.set_visible(style == "block")
            sp.set_linewidth(0.25 * S.g); sp.set_color("black")
    return top


def aligned_blocks(ctx):
    """(read_start, read_end, ref_start) for every aligned block of the primary
    alignment, in the coordinates the quad plot uses for each block."""
    a = ctx.primary
    out, q, r = [], 0, a.reference_start
    for op, ln in (a.cigartuples or []):
        if op in (0, 7, 8):
            qs = _seq_x(q, a, ctx)
            out.append((qs, qs + ln, r - ctx.win_start))
            q += ln
            r += ln
        elif op == 1:
            q += ln
        elif op in (2, 3):
            r += ln
        elif op in (4, 5):
            q += ln
    return out


def ref_at(blocks, q):
    """The reference position the aligner gives read position `q`.

    Inside an aligned block it is that block's own mapping. Inside an INSERTION it is the
    insertion site -- the reference position the aligner had reached when it took those
    bases from the read."""
    if not blocks:
        return 0
    prev = blocks[0][2]
    for qs, qe, rs in blocks:
        if q < qs:
            return prev
        if q < qe:
            return rs + (q - qs)
        prev = rs + (qe - qs)
    return prev


def insertion_sites(ctx, read_intervals):
    """Where each annotated read interval sits on the reference, when it sits at a POINT.

    A read interval that lies wholly inside an insertion has no reference span: the
    aligner took those bases without advancing along the reference. That single position
    is the insertion site, and it is the only place the reference block can show an
    insertion at all. Which of an INS row's two read intervals is the inserted copy and
    which is the donor's copy is not in the annotation columns -- the alignment says so,
    by projecting one of them to zero width."""
    blocks = aligned_blocks(ctx)
    R = len(ctx.ref_seq)
    out = []
    for a, b in read_intervals:
        x0, x1 = ref_at(blocks, a - R), ref_at(blocks, b - R)
        if x1 - x0 < 1 and 0 <= x0 <= R:
            out.append(x0)
    return sorted(set(out))


def _annotation(ax, ctx, lines, S):
    """Black rectangles around the annotated intervals, in every quadrant that holds one.

    An annotated interval -- the donor region, the inserted segment, the deleted block --
    is a stretch of sequence, and what it produces in a dot plot is a diagonal. The box is
    what that diagonal is FOR.

    Only the two self-comparison quadrants are boxed: reference x reference and read x
    read, each on its own diagonal. Whatever a cross quadrant is given -- a rectangle per
    (reference, read) interval pair, or each read interval at the reference span the
    aligner gives it -- says something the two diagonal boxes and the data between them
    already say, at the cost of two more marks per interval.

    An event that is a JUNCTION rather than a segment gets a dotted cross-hair, drawn
    inside the block it belongs to. There are two, and they are each other's mirror: a
    deletion is a junction in the read (the reference block boxes what is missing, the
    read has only a step in a diagonal to show for it), and an insertion is a junction in
    the reference (the read block boxes the inserted copy, the reference never held it).

    Returns (boxes, cross-hairs)."""
    R, Q = len(ctx.ref_seq), len(ctx.read_seq)
    ref, read = lines.intervals(ctx)
    rects = [(a, b, a, b) for a, b in ref] + [(a, b, a, b) for a, b in read]
    for x0, x1, y0, y1 in rects:
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False,
                               edgecolor=COL_BOX, linewidth=BOX_LW * S.g, zorder=6))

    drawn = {p for a, b in read for p in (a, b)}
    dels = [p for p in sorted({R + q for q in lines.read})
            if p not in drawn and R <= p <= R + Q]
    ins = insertion_sites(ctx, read)
    for marks, lo, hi in ((dels, R, R + Q), (ins, 0, R)):
        for m in marks:
            ax.vlines(m, lo, hi, colors=COL_BOX, lw=BOX_LW * S.g, ls=":", zorder=6)
            ax.hlines(m, lo, hi, colors=COL_BOX, lw=BOX_LW * S.g, ls=":", zorder=6)
    return len(rects), (len(ins), len(dels))


def block_labels(ctx, coords=False):
    """What the two halves of each axis are called: the reference window by its
    coordinates, the read by its length. `coords` keeps the unit of the tick labels,
    which only the coordinate decoration draws."""
    ref = f"{ctx.chrom}:{ctx.win_start + 1:,}-{ctx.win_end:,}"
    read = f"Read ({len(ctx.read_seq):,} bp)"
    return (f"{ref}  (Mb)", "Read (kb)") if coords else (ref, read)


def _quad_figure(ctx, params, lines, fwd, rev, n, track, S):
    """One drawing of the quad plot, at scale `S`. Returns (fig, strip axes, counts)."""
    k = params.kmer
    R, Q = len(ctx.ref_seq), len(ctx.read_seq)

    mpl.rcParams.update(S.rc())
    box = S.panel_mm * MM
    style = params.monomer_style if params.monomer_style == "block" else "lollipop"
    gap = STRIP_GAP_MM * MM * S.g
    if style == "block" or track is None:
        strip = STRIP_BLOCK_MM * MM * S.g
    else:   # exactly what the marks need: the stick, plus the circle sitting on its end
        strip = (STEM_MM * MM * S.g + _dot_size(S, len(track.units)) / 72 / 2
                 + 0.05 * MM * S.g)
    if track is None:      # inches of margin, room for ticks; they hold text, so they
        ml, mr, mb, mt = (0.46 * S.t, 0.08 * S.t, 0.46 * S.t, 0.44 * S.t)   # follow it
    else:                                          # no tick labels, but strips instead
        ml, mr = 0.17 * S.t, 0.08 * S.t + strip + gap
        mb, mt = 0.17 * S.t, 0.44 * S.t + strip + gap
    fw, fh = box + ml + mr, box + mb + mt
    fig = plt.figure(figsize=(fw, fh))
    ax = fig.add_axes([ml / fw, mb / fh, box / fw, box / fh])

    if lines and params.annot_style in ("lines", "both"):
        for v in lines.marks(ctx):
            ax.axvline(v, color=COL_LINE, lw=0.25 * S.g, ls=":", zorder=1)
            ax.axhline(v, color=COL_LINE, lw=0.25 * S.g, ls=":", zorder=1)

    bands = identity_bands(ctx, R)
    a0, b0, a1, b1 = fwd
    if a0.size:
        xy = np.stack([np.column_stack([a0, b0]),
                       np.column_stack([a1 + k - 1, b1 + k - 1])], axis=1)
        run_d, run_hi = b0 - a0, a1 + k
        need = np.maximum(1, (run_hi - a0) * MIN_OVERLAP)
        is_main = np.zeros(a0.size, dtype=bool)
        for d, lo, hi in extend_bands(bands, a0, run_hi, run_d):
            overlap = np.minimum(run_hi, hi) - np.maximum(a0, lo)
            is_main |= (np.abs(run_d - d) <= DIAG_TOL) & (overlap >= need)
        for sel, col, z in ((~is_main, params.colour_ext, 2),
                            (is_main, params.colour_main, 3)):
            if sel.any():
                ax.add_collection(LineCollection(xy[sel], colors=col,
                                                 linewidths=0.25 * S.g, zorder=z,
                                                 rasterized=int(sel.sum()) > 20000))
    a0, b0, a1, b1 = rev
    if a0.size:
        xy = np.stack([np.column_stack([a0, b0 + k - 1]),
                       np.column_stack([a1 + k - 1, b1])], axis=1)
        ax.add_collection(LineCollection(xy, colors=COL_REV, linewidths=0.25 * S.g,
                                         zorder=4, rasterized=xy.shape[0] > 20000))

    n_box, (n_ins, n_del) = 0, (0, 0)
    if lines and params.annot_style in ("box", "both"):
        n_box, (n_ins, n_del) = _annotation(ax, ctx, lines, S)

    ax.set_xlim(0, n); ax.set_ylim(0, n)
    ax.axvline(R, color="black", lw=0.5 * S.g, zorder=5)
    ax.axhline(R, color="black", lw=0.5 * S.g, zorder=5)
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_linewidth(0.5 * S.g); sp.set_color("black")
    ax.tick_params(width=0.5 * S.g, color="black")

    strip_ax = None
    if track is not None:
        ax.set_xticks([]); ax.set_yticks([])
        strip_ax = _strips(fig, (ml, mb, box, gap, strip, fw, fh), R, n, track, style, S)
    else:
        _coordinate_ticks(fig, ax, ctx, R, Q, S)

    off_x, off_y = ((-6, -8) if track is not None else (-17, -21))
    for pos, lab in zip((R / 2, R + Q / 2), block_labels(ctx, coords=track is None)):
        ax.annotate(lab, xy=(pos, 0), xycoords=("data", "axes fraction"),
                    xytext=(0, off_x * S.t), textcoords="offset points",
                    ha="center", va="top", fontsize=AXLAB_PT * S.t, color=COL_AXLAB)
        ax.annotate(lab, xy=(0, pos), xycoords=("axes fraction", "data"),
                    xytext=(off_y * S.t, 0), textcoords="offset points", rotation=90,
                    ha="right", va="center", fontsize=AXLAB_PT * S.t, color=COL_AXLAB)

    title = (f"{ctx.read_id}\n{ctx.window}  (strand {ctx.strand})\n"
             f"ref {R:,} + read {Q:,} bp  ·  k={k}, min_seg={params.min_seg}")
    if track is not None:                          # its own line: the title sets the
        title += (f"\n{track.n_full} monomers of {track.period} bp, "   # figure width
                  f"{track.n_groups} groups at {int(params.monomer_cut * 100)}% identity"
                  f"\nphase: {track.phase}"
                  + (f"  ·  {track.n_nonsatellite} non-satellite"
                     if track.n_nonsatellite else ""))
    if n_box:
        title += "\nbox: annotated donor, inserted or deleted segment"
    notes = (["the insertion site in the reference"] if n_ins else []) + \
            (["the deletion junction in the read"] if n_del else [])
    if notes:
        title += "\ndotted: " + " and ".join(notes)
    (strip_ax or ax).set_title(title, fontsize=BASE_PT * S.t, linespacing=1.6)
    return fig, (n_box, n_ins, n_del)


def quad(ctx, params, out_stem, lines=None, formats=("png", "pdf")):
    """Draw the quad plot. Returns (paths, stats).

    Each format is drawn at its own size: the PDF at `--panel-mm` for placing in a figure,
    the PNG larger and with larger text for looking at on a screen. The k-mer comparison
    and the monomer annotation are done once and shared."""
    R, Q = len(ctx.ref_seq), len(ctx.read_seq)
    fwd, rev, n = _self_compare(ctx, params)
    track = mono.annotate(ctx, period=params.monomer_period, cut=params.monomer_cut,
                          consensus=params.monomer_consensus) if params.monomer else None

    paths, counts = [], (0, 0, 0)
    for fmt in formats:
        S = sizes_for(fmt, params)
        fig, counts = _quad_figure(ctx, params, lines, fwd, rev, n, track, S)
        path = f"{out_stem}.quad.{fmt}"
        fig.savefig(path, format=fmt, dpi=params.dpi)
        plt.close(fig)
        paths.append(path)
    return paths, dict(n_fwd=int(fwd[0].size), n_rev=int(rev[0].size), total_bp=n,
                       monomer=track, n_boxes=counts[0], n_marks=counts[1] + counts[2])


def _coordinate_ticks(fig, ax, ctx, R, Q, S):
    """The default decoration: genomic Mb on the reference block, kb on the read."""
    ticks, labels, n_ref = _ticks(ctx, R, Q)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(mpl.ticker.FixedLocator(ticks))
        axis.set_major_formatter(mpl.ticker.FixedFormatter(labels))
    # The reference block's last label and the read block's first label sit on the same
    # coordinate (the divider); align them outward so they meet rather than overprint.
    ends = {0: "start", n_ref - 1: "end", n_ref: "start", len(ticks) - 1: "end"}
    pad_r = ScaledTranslation(0.022 * S.g, 0, fig.dpi_scale_trans)   # ~0.5 mm
    pad_l = ScaledTranslation(-0.022 * S.g, 0, fig.dpi_scale_trans)
    for i, t in enumerate(ax.get_xticklabels()):
        if i in ends:
            start = ends[i] == "start"
            t.set_horizontalalignment("left" if start else "right")
            t.set_transform(t.get_transform() + (pad_r if start else pad_l))
    for i, t in enumerate(ax.get_yticklabels()):
        if i in ends:
            t.set_verticalalignment("bottom" if ends[i] == "start" else "top")
