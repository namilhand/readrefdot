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
from matplotlib.transforms import ScaledTranslation

from .kmer import (build_index, filter_min_length, kmer_codes, match,
                   merge_segments, revcomp)

MM = 1 / 25.4
PANEL_MM = 45.0                 # default plot box, excluding title and labels
TICK_STEP = 5000

COL_MAIN = "#000000"            # diagonals the aligner placed the read on
COL_EXT = "#B2B2B2"             # grey70 - every other diagonal
COL_REV = "#D55E00"             # vermillion - reverse-complement matches
COL_LINE = "#0072B2"            # blue - annotation guide lines
COL_LAB = "#444444"
DIAG_TOL = 3                    # bp slack when matching a run to an alignment diagonal
MIN_BLOCK = 20                  # smallest aligned block that contributes an identity band
MIN_OVERLAP = 0.5               # a run must lie this far inside a band to count as main
EXTEND_GAP = 200                # bp: bands grow along their diagonal through runs this close

STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
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


@dataclass
class Params:
    kmer: int = 20
    min_seg: int = 170          # ~one CEN178 satellite unit (units vary 177-179 bp)
    merge_gap: int = None       # default k+1: bridges a single substitution
    colour_main: str = COL_MAIN
    colour_ext: str = COL_EXT
    panel_mm: float = PANEL_MM

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


def quad(ctx, params, out_stem, lines=None, formats=("png", "pdf")):
    """Draw the quad plot. Returns (paths, stats)."""
    k = params.kmer
    R, Q = len(ctx.ref_seq), len(ctx.read_seq)
    fwd, rev, n = _self_compare(ctx, params)

    mpl.rcParams.update(STYLE)
    box = params.panel_mm * MM
    ml, mr, mb, mt = 0.46, 0.08, 0.46, 0.44        # inches of margin
    fw, fh = box + ml + mr, box + mb + mt
    fig = plt.figure(figsize=(fw, fh))
    ax = fig.add_axes([ml / fw, mb / fh, box / fw, box / fh])

    if lines:
        for v in lines.marks(ctx):
            ax.axvline(v, color=COL_LINE, lw=0.25, ls=":", zorder=1)
            ax.axhline(v, color=COL_LINE, lw=0.25, ls=":", zorder=1)

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
                ax.add_collection(LineCollection(xy[sel], colors=col, linewidths=0.25,
                                                 zorder=z,
                                                 rasterized=int(sel.sum()) > 20000))
    a0, b0, a1, b1 = rev
    if a0.size:
        xy = np.stack([np.column_stack([a0, b0 + k - 1]),
                       np.column_stack([a1 + k - 1, b1])], axis=1)
        ax.add_collection(LineCollection(xy, colors=COL_REV, linewidths=0.25, zorder=4,
                                         rasterized=xy.shape[0] > 20000))

    ax.set_xlim(0, n); ax.set_ylim(0, n)
    ax.axvline(R, color="black", lw=0.5, zorder=5)
    ax.axhline(R, color="black", lw=0.5, zorder=5)
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_linewidth(0.5); sp.set_color("black")
    ax.tick_params(width=0.5, color="black")

    ticks, labels, n_ref = _ticks(ctx, R, Q)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(mpl.ticker.FixedLocator(ticks))
        axis.set_major_formatter(mpl.ticker.FixedFormatter(labels))
    # The reference block's last label and the read block's first label sit on the same
    # coordinate (the divider); align them outward so they meet rather than overprint.
    ends = {0: "start", n_ref - 1: "end", n_ref: "start", len(ticks) - 1: "end"}
    pad_r = ScaledTranslation(0.022, 0, fig.dpi_scale_trans)     # ~0.5 mm
    pad_l = ScaledTranslation(-0.022, 0, fig.dpi_scale_trans)
    for i, t in enumerate(ax.get_xticklabels()):
        if i in ends:
            start = ends[i] == "start"
            t.set_horizontalalignment("left" if start else "right")
            t.set_transform(t.get_transform() + (pad_r if start else pad_l))
    for i, t in enumerate(ax.get_yticklabels()):
        if i in ends:
            t.set_verticalalignment("bottom" if ends[i] == "start" else "top")

    for pos, lab in ((R / 2, f"{ctx.chrom} (Mb)"), (R + Q / 2, "read (kb)")):
        ax.annotate(lab, xy=(pos, 0), xycoords=("data", "axes fraction"),
                    xytext=(0, -17), textcoords="offset points",
                    ha="center", va="top", fontsize=5.5, color=COL_LAB)
        ax.annotate(lab, xy=(0, pos), xycoords=("axes fraction", "data"),
                    xytext=(-21, 0), textcoords="offset points", rotation=90,
                    ha="right", va="center", fontsize=5.5, color=COL_LAB)

    ax.set_title(f"{ctx.read_id}\n{ctx.window}  (strand {ctx.strand})\n"
                 f"ref {R:,} + read {Q:,} bp  ·  k={k}, min_seg={params.min_seg}",
                 fontsize=5, linespacing=1.6)

    paths = []
    for fmt in formats:
        path = f"{out_stem}.quad.{fmt}"
        fig.savefig(path, format=fmt)
        paths.append(path)
    plt.close(fig)
    return paths, dict(n_fwd=int(fwd[0].size), n_rev=int(rev[0].size), total_bp=n)
