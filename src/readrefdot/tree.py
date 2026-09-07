"""Neighbour-joining tree of the monomers in one plot.

The dot plot says where the satellite units are; the tree says how they are related.
Tips are the units of that plot, coloured by the same similarity group as the axis
annotation, so a colour in the tree and a colour on the axis are the same variant.

The distance is alignment-free (1 - the shared-k-mer identity estimate used for
grouping), not a substitution model, so the scale bar is labelled as distance rather
than substitutions per site.
"""

import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

from .monomer import PALETTE, UNSET
from .plot import MM, STYLE

TIP_MIN, TIP_MAX = 0.8, 2.5        # pt, marker diameter
SCALE_STEPS = (0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2)


def neighbour_joining(D):
    """Standard NJ. Returns (children, root) with children[node] = [(child, length)].

    Nodes 0..n-1 are the tips, in the order of D."""
    n = D.shape[0]
    if n < 3:
        raise ValueError("neighbour joining needs at least 3 tips")
    size = 2 * n - 1
    d = np.zeros((size, size))
    d[:n, :n] = D
    active = list(range(n))
    children = {}
    nxt = n
    while len(active) > 2:
        idx = np.array(active)
        m = idx.size
        sub = d[np.ix_(idx, idx)]
        r = sub.sum(axis=1)
        q = (m - 2) * sub - r[:, None] - r[None, :]
        np.fill_diagonal(q, np.inf)
        i, j = np.unravel_index(np.argmin(q), q.shape)
        li = sub[i, j] / 2 + (r[i] - r[j]) / (2 * (m - 2))
        lj = sub[i, j] - li
        u = nxt
        nxt += 1
        children[u] = [(int(idx[i]), max(0.0, float(li))),
                       (int(idx[j]), max(0.0, float(lj)))]
        new = (sub[i, :] + sub[j, :] - sub[i, j]) / 2
        d[u, idx] = new
        d[idx, u] = new
        d[u, u] = 0.0
        active = [a for a in active if a not in (idx[i], idx[j])] + [u]
    a, b = active
    half = float(d[a, b]) / 2
    children[nxt] = [(a, half), (b, half)]     # midpoint root: only a display choice
    return children, nxt


def _leaf_counts(children, root):
    counts, order = {}, []
    stack = [(root, False)]
    while stack:
        node, done = stack.pop()
        if done:
            counts[node] = sum(counts[c] for c, _ in children[node])
            continue
        if node not in children:
            counts[node] = 1
            order.append(node)
            continue
        stack.append((node, True))
        for c, _ in children[node]:
            stack.append((c, False))
    return counts


def equal_angle_layout(children, root):
    """Felsenstein's equal-angle layout: each subtree gets a wedge proportional to how
    many tips it holds, and a branch is drawn at its wedge's mid-angle with a length
    equal to the branch length. Unrooted in appearance, which is how a tree of variants
    with no outgroup should be read."""
    counts = _leaf_counts(children, root)
    pos = {root: (0.0, 0.0)}
    stack = [(root, 0.0, 2 * np.pi)]
    while stack:
        node, a0, a1 = stack.pop()
        kids = children.get(node)
        if not kids:
            continue
        total = sum(counts[c] for c, _ in kids)
        start = a0
        for c, blen in kids:
            wedge = (a1 - a0) * counts[c] / total
            mid = start + wedge / 2
            x, y = pos[node]
            pos[c] = (x + blen * np.cos(mid), y + blen * np.sin(mid))
            stack.append((c, start, start + wedge))
            start += wedge
    return pos


def _scale_value(span):
    """A round scale-bar length, about a quarter of the tree's width."""
    target = span / 4 if span > 0 else SCALE_STEPS[0]
    return min(SCALE_STEPS, key=lambda s: abs(np.log(s / target)) if s > 0 else np.inf)


def draw(track, out_stem, panel_mm=45.0, formats=("png", "pdf"), title=None):
    """Write <out_stem>.tree.png/.pdf. Returns the paths, or [] if there is no tree."""
    ident = track.identity
    if ident is None or ident.shape[0] < 3:
        return []
    full = [u for u in track.units if not u.partial]
    dist = np.clip(1.0 - ident, 0.0, None)
    np.fill_diagonal(dist, 0.0)
    children, root = neighbour_joining(dist)
    pos = equal_angle_layout(children, root)

    segs = [[pos[node], pos[c]] for node, kids in children.items() for c, _ in kids]
    xy = np.array([pos[i] for i in range(len(full))])
    cols = [UNSET if u.group < 0 or u.group >= len(PALETTE) else PALETTE[u.group]
            for u in full]

    mpl.rcParams.update(STYLE)
    box = panel_mm * MM
    ml, mr, mb, mt = 0.06, 0.06, 0.22, 0.40 if title else 0.06
    fw, fh = box + ml + mr, box + mb + mt
    fig = plt.figure(figsize=(fw, fh))
    ax = fig.add_axes([ml / fw, mb / fh, box / fw, box / fh])
    ax.add_collection(LineCollection(segs, colors="black", linewidths=0.25, zorder=1))

    # Frame on the drawing's bounding box, not on the mean of the node positions: an
    # equal-angle layout is not centred on its root, and centring on the mean clips the
    # longest spokes.
    allxy = np.array(list(pos.values()))
    lo, hi = allxy.min(axis=0), allxy.max(axis=0)
    span = float(max(hi - lo)) or 1.0
    cx, cy = (lo + hi) / 2
    # Tips spread over two dimensions here, so they can be larger than the axis dots.
    tip = float(np.clip(panel_mm / np.sqrt(max(len(full), 1)) * 0.25 / 25.4 * 72,
                        TIP_MIN, TIP_MAX))
    ax.scatter(xy[:, 0], xy[:, 1], s=tip ** 2, c=cols, linewidths=0, zorder=2)

    pad = span * 0.06
    width = span + 2 * pad
    ax.set_xlim(cx - width / 2, cx + width / 2)
    ax.set_ylim(cy - width / 2, cy + width / 2)
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    bar = _scale_value(span)
    frac = bar / width
    x0, y0 = 0.03, -0.04                         # axes fraction, just below the frame
    ax.plot([x0, x0 + frac], [y0, y0], transform=ax.transAxes,
            color="black", lw=0.5, clip_on=False)
    ax.text(x0 + frac / 2, y0 - 0.03, f"{bar:g} distance",
            transform=ax.transAxes, ha="center", va="top", fontsize=5, color="#444444")
    if title:
        ax.set_title(title, fontsize=5, linespacing=1.6)

    paths = []
    for fmt in formats:
        path = f"{out_stem}.tree.{fmt}"
        fig.savefig(path, format=fmt)
        paths.append(path)
    plt.close(fig)
    return paths
