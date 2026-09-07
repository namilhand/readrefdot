"""satdivplot - pairwise divergence between the satellite monomers of one sequence.

The quad dot plot shows where sequence recurs; this shows how far apart the copies are.
Each panel is a self-comparison: the array is cut into CEN178 monomers (the same tiling
and the same grouping readrefdot uses), every monomer is compared with every other, and
the resulting n x n matrix of percent divergence is drawn as a heat map, monomers in
array order on both axes.

Higher-order repeat structure is what this makes visible. In an array built of a
repeating cassette of m monomers, monomer i and monomer i+m are near-identical while
their neighbours are not, so the matrix carries a chequer of low-divergence cells at
multiples of the HOR period -- off-diagonals the dot plot only hints at.

One page holds two panels: the reference window on the left, the read on the right. Both
are tiled and grouped together, so a colour on the top/right strips means the same
monomer family in both.

Divergence is measured in CONSENSUS COLUMNS. Each monomer is aligned to the consensus
and projected onto its 178 columns, so an indel inside one monomer shifts nothing
downstream and two monomers are always compared base-for-corresponding-base.
"""

from dataclasses import dataclass

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from . import monomer as mono
from .plot import MM, STEM_MM, STYLE, _dot_size, _lollipops

PANEL_MM = 25.0            # default heat map box, per block
GAP_MM = 0.4               # panel to strip
INTER_MM = 5.0             # between the reference group and the read group
CB_GAP_MM = 3.5            # strip to colour bar
CB_W_MM = 1.6
CHUNK = 64                 # rows of the matrix computed at once
DPI = 600                  # the raster output: 25 mm holds ~120 cells, so 300 is coarse
BOX_LW = 0.5               # the black square drawn around an annotated interval
CMAP = "RdYlBu_r"          # diverging: blue = alike, red = far apart
LO_PCT, HI_PCT = 2.0, 98.0  # percentiles of the off-diagonal that set the colour range
COL_LAB = "#444444"


@dataclass
class Params:
    panel_mm: float = PANEL_MM
    cmap: str = CMAP
    vmin: float = None         # % divergence at the blue end
    vmax: float = None         # % divergence at the red end
    center: float = None       # % divergence the map is neutral at (default: the median)
    dpi: int = DPI             # raster resolution; the PDF stays vector either way
    monomer_period: int = None
    monomer_cut: float = mono.DEFAULT_CUT
    monomer_consensus: str = mono.CEN178


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


def block_matrix(S, track, block, cons):
    """(units, colours, divergence matrix) for one block, in array order.

    Only whole satellite monomers are compared -- the same ones the grouping used. A
    partial unit at the edge of the window is a fragment of a monomer, and a
    non-satellite stretch is not a monomer at all; either would be a spurious row."""
    units = [u for u in track.units
             if u.block == block and not u.partial and u.satellite]
    P = []
    keep = []
    for u in units:
        p = project(S[u.start:u.end], cons)
        if p is not None:
            P.append(p)
            keep.append(u)
    if not P:
        return keep, [], np.zeros((0, 0))
    return keep, [_colour(u) for u in keep], divergence(np.array(P, dtype=np.uint8))


def _index_at(units, pos):
    """Fractional monomer index of a base position: unit 3 half way through is 2.5.

    Interpolating inside the unit rather than snapping to it puts the edge of a drawn box
    on the base the annotation names, not on the nearest monomer boundary."""
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


def _intervals(ctx, lines, block):
    """The annotated intervals for one block, in the coordinates the units are tiled in.

    `ref-lines`/`read-lines` hold the same intervals readrefdot draws: for an INS the
    donor region on the reference, and in the read both the donor's copy and the inserted
    segment; for a DEL the deleted block on the reference. `Lines.intervals` pairs the
    endpoints (see `annotate.pair_up`) and clips each to its own block."""
    if not lines:
        return []
    ref, read = lines.intervals(ctx)
    return ref if block == "ref" else read


def _draw_boxes(ax, units, spans):
    """A black square around each annotated interval, on the diagonal.

    A run of monomers that recur elsewhere in the array draws a diagonal line; the square
    is what that line is FOR -- the donor segment, or the inserted copy of it."""
    n = 0
    for a, b in spans:
        if b <= units[0].start or a >= units[-1].end:
            continue                    # the interval is outside this block's monomers
        i, j = _index_at(units, a), _index_at(units, b)
        if j - i < 0.5:                 # a point, not a segment (a deletion junction)
            continue
        ax.add_patch(Rectangle((i, i), j - i, j - i, fill=False, edgecolor="black",
                               linewidth=BOX_LW, zorder=5, clip_on=True))
        n += 1
    return n


def _limits(mats, p):
    """One colour range for both panels, so the two are directly comparable.

    The scale is diverging and is centred on the MEDIAN divergence of the array, not on
    the middle of the range: the neutral colour then means "as different as two monomers
    of this array typically are", blue means closer kin than that and red more distant.
    An array whose monomers all sit near 5% and one that spans 1-9% therefore read the
    same way, which is what makes a repeating pattern of blue cells legible as HOR
    structure rather than as the array's overall age."""
    parts = [m[~np.eye(m.shape[0], dtype=bool)].ravel() for m in mats if m.size]
    if not parts:
        return 0.0, 1.0
    off = np.concatenate(parts)
    mid = float(p.center) if p.center is not None else float(np.median(off))
    lo = float(p.vmin) if p.vmin is not None else float(np.percentile(off, LO_PCT))
    hi = float(p.vmax) if p.vmax is not None else float(np.percentile(off, HI_PCT))
    if p.vmin is None and p.vmax is None:      # symmetric about the centre, but never
        r = max(mid - lo, hi - mid)            # reaching below zero divergence
        lo = max(0.0, mid - r)
        hi = mid + (mid - lo)
    return (lo, hi) if hi > lo else (lo, lo + 1.0)


def _panel(fig, geom, D, cols, cmap, lo, hi, label, span, dot, stem,
           units=None, spans=()):
    """One heat map with its monomer lollipops along the top and the right.

    The axes run left to right and BOTTOM TO TOP, so monomer 1 is at the bottom-left
    corner and the array reads outwards in both directions.

    `span` is the number of monomer slots the box holds, and is the same for both panels
    so a cell is the same size in each and the two can be compared directly. A block with
    fewer monomers than `span` simply leaves the far end of the box empty -- the room an
    insertion takes up in the other block."""
    x0, y0, box, gap, strip, fw, fh = geom
    n = D.shape[0]
    ax = fig.add_axes([x0 / fw, y0 / fh, box / fw, box / fh])
    im = None
    if n:
        im = ax.imshow(D, cmap=cmap, vmin=lo, vmax=hi, extent=(0, n, 0, n),
                       origin="lower", interpolation="nearest", aspect="auto")
    n_box = _draw_boxes(ax, units, spans) if (n and spans) else 0
    ax.set_xlim(0, span); ax.set_ylim(0, span)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_linewidth(0.5); sp.set_color("black")

    centres = np.arange(n) + 0.5
    top = fig.add_axes([x0 / fw, (y0 + box + gap) / fh, box / fw, strip / fh])
    right = fig.add_axes([(x0 + box + gap) / fw, y0 / fh, strip / fw, box / fh])
    _lollipops(top, centres, cols, dot, True, stem)
    _lollipops(right, centres, cols, dot, False, stem)
    top.set_xlim(0, span); top.set_ylim(0, 1)
    right.set_xlim(0, 1); right.set_ylim(0, span)   # matches the heat map's y direction
    for a in (top, right):
        a.set_xticks([]); a.set_yticks([])
        for name, sp in a.spines.items():
            sp.set_visible(name in (("bottom",) if a is top else ("left",)))
            sp.set_linewidth(0.25); sp.set_color("black")

    ax.annotate(label, xy=(0.5, 0), xycoords="axes fraction", xytext=(0, -5),
                textcoords="offset points", ha="center", va="top", fontsize=5.5,
                color=COL_LAB)
    return ax, im, n_box


def draw(ctx, params, out_stem, lines=None, formats=("pdf", "png")):
    """Draw the two-panel divergence figure. Returns (paths, stats)."""
    S = ctx.ref_seq + ctx.read_seq
    track = mono.annotate(ctx, period=params.monomer_period, cut=params.monomer_cut,
                          consensus=params.monomer_consensus)
    if track is None:
        raise ValueError("not a tandem satellite array: no monomers to compare")
    cons = params.monomer_consensus or S[track.units[0].start:track.units[0].end]
    if params.monomer_consensus:
        # annotate() may have phased the array on the reverse strand; project onto the
        # same orientation it tiled with, or every unit would align back to front.
        cons = (mono.revcomp(params.monomer_consensus) if "reverse" in track.phase
                else params.monomer_consensus)

    panels = [(name, *block_matrix(S, track, name, cons)) for name in ("ref", "read")]
    lo, hi = _limits([p[3] for p in panels], params)

    mpl.rcParams.update(STYLE)
    box = params.panel_mm * MM
    # Both boxes hold the same number of monomer slots, so a cell is the same size in
    # each: the shorter block leaves the far end of its box empty rather than stretching
    # to fill it, and the gap is exactly the length the other block has gained.
    span = max(p[3].shape[0] for p in panels) or 1
    dot = _dot_size(params.panel_mm, span)
    stem_mm = STEM_MM * MM
    strip = stem_mm + dot / 72 / 2 + 0.05 * MM      # the stick plus half a circle
    gap, inter = GAP_MM * MM, INTER_MM * MM
    cb_gap, cb_w = CB_GAP_MM * MM, CB_W_MM * MM
    ml, mb, mt, mr = 0.10, 0.16, 0.50, 0.34
    group = box + gap + strip
    fw = ml + 2 * group + inter + cb_gap + cb_w + mr
    fh = mb + group + mt
    fig = plt.figure(figsize=(fw, fh))

    labels = {"ref": f"reference  {ctx.chrom}", "read": "read"}
    im, n_box = None, 0
    for i, (name, units, cols, D) in enumerate(panels):
        x0 = ml + i * (group + inter)
        _, got, nb = _panel(fig, (x0, mb, box, gap, strip, fw, fh), D, cols, params.cmap,
                            lo, hi, f"{labels[name]}  ({D.shape[0]} monomers)", span,
                            dot, stem_mm / strip, units,
                            _intervals(ctx, lines, name))
        im = im or got
        n_box += nb

    cax = fig.add_axes([(ml + 2 * group + inter + cb_gap) / fw, mb / fh, cb_w / fw,
                        box / fh])
    cb = fig.colorbar(im, cax=cax)
    cb.outline.set_linewidth(0.4)
    cax.tick_params(width=0.4, length=1.6, labelsize=4.5, pad=1.2)
    cax.set_ylabel("monomer-pair divergence (%)", fontsize=5, color=COL_LAB, labelpad=2)

    n_ref, n_read = panels[0][3].shape[0], panels[1][3].shape[0]
    title = (f"{ctx.read_id}\n{ctx.window}  (strand {ctx.strand})\n"
             f"{n_ref} + {n_read} monomers of {track.period} bp, "
             f"{track.n_groups} groups at {int(params.monomer_cut * 100)}% identity"
             f"\nphase: {track.phase}"
             + (f"  ·  {track.n_nonsatellite} non-satellite" if track.n_nonsatellite
                else "")
             + ("\nblack box: annotated donor, inserted or deleted segment" if n_box else ""))
    fig.text(ml / fw, 1 - 0.012, title, ha="left", va="top", fontsize=5,
             linespacing=1.6)

    paths = []
    for fmt in formats:
        path = f"{out_stem}.satdiv.{fmt}"
        fig.savefig(path, format=fmt, dpi=params.dpi)
        paths.append(path)
    plt.close(fig)
    return paths, dict(n_ref=n_ref, n_read=n_read, vmin=lo, vmax=hi, monomer=track,
                       n_boxes=n_box,
                       matrices={n: D for n, _, _, D in panels})


def write_tsv(ctx, track, name, units, D, path):
    """The matrix itself: one row per monomer, one column per monomer."""
    with open(path, "w") as fh:
        head = "\t".join(f"{i + 1}" for i in range(len(units)))
        fh.write(f"block\tindex\tstart\tend\tgroup\t{head}\n")
        R = len(ctx.ref_seq)
        for i, u in enumerate(units):
            b = u.start if name == "ref" else u.start - R
            vals = "\t".join(f"{v:.3f}" for v in D[i])
            fh.write(f"{name}\t{i + 1}\t{b}\t{b + u.length}\t"
                     f"{u.group if u.group >= 0 else '.'}\t{vals}\n")
    return path
