"""The monomer dendrogram — and, optionally, a neighbour-joining tree.

The dendrogram is the SAME clustering that colours the axis: one average-linkage tree,
cut once. So a colour is a branch, the cut is a line you can see, and the picture cannot
disagree with the annotation. Neighbour joining is kept for when a rate-corrected,
unrooted topology is wanted (`--tree-method nj`); it is a second clustering of the same
distances and will not agree with the colours exactly.

Neighbour-joining tree of the monomers in one plot.

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


def _dendrogram_layout(merges, n):
    """Leaf order and node coordinates for a rectangular dendrogram.

    Leaves are ordered by an in-order walk of the merge tree, which is what keeps the
    branches from crossing."""
    kids = {}
    nxt = n
    for a, b, h in merges:
        kids[nxt] = (a, b, h)
        nxt += 1
    root = nxt - 1
    order, stack = [], [root]
    while stack:                                  # iterative in-order walk
        node = stack.pop()
        if node < n:
            order.append(node)
        else:
            a, b, _ = kids[node]
            stack.extend((b, a))                  # a first
    xpos = {leaf: i for i, leaf in enumerate(order)}
    height = {leaf: 0.0 for leaf in range(n)}
    for node in range(n, nxt):
        a, b, h = kids[node]
        xpos[node] = (xpos[a] + xpos[b]) / 2
        height[node] = h
    return order, xpos, height, kids, root


def draw_dendrogram(track, out_stem, panel_mm=45.0, formats=("png", "pdf"), title=None):
    """Write <out_stem>.dendrogram.png/.pdf: the grouping tree, cut line and all."""
    n = 0 if track.identity is None else track.identity.shape[0]
    if n < 3 or not track.merges:
        return []
    full = [u for u in track.units if not u.partial and u.satellite]
    order, xpos, height, kids, root = _dendrogram_layout(track.merges, n)

    mpl.rcParams.update(STYLE)
    box = panel_mm * MM
    ml, mr, mb, mt = 0.30, 0.06, 0.20, 0.44 if title else 0.06
    fw, fh = box + ml + mr, box + mb + mt
    fig = plt.figure(figsize=(fw, fh))
    ax = fig.add_axes([ml / fw, mb / fh, box / fw, box / fh])

    def colour_of(node):
        """A branch takes the group's colour while it is still inside that group."""
        if height[node] > track.cut:
            return "#000000"
        leaf = node
        while leaf >= n:
            leaf = kids[leaf][0]
        g = full[leaf].group
        return PALETTE[g] if 0 <= g < len(PALETTE) else UNSET

    segs, cols = [], []
    for node in range(n, n + len(track.merges)):
        a, b, h = kids[node]
        for child in (a, b):
            segs.append([(xpos[child], height[child]), (xpos[child], h)])
            cols.append(colour_of(child))
        segs.append([(xpos[a], h), (xpos[b], h)])
        cols.append(colour_of(node))
    ax.add_collection(LineCollection(segs, colors=cols, linewidths=0.35, zorder=2))

    tip = float(np.clip(panel_mm / max(n, 1) * 0.9 / 25.4 * 72, TIP_MIN, TIP_MAX))
    tcol = [PALETTE[full[i].group] if 0 <= full[i].group < len(PALETTE) else UNSET
            for i in range(n)]
    ax.scatter([xpos[i] for i in range(n)], np.zeros(n), s=tip ** 2,
               c=tcol, linewidths=0, zorder=3, clip_on=False)

    top = max(height.values()) or 1.0
    ax.axhline(track.cut, color="#444444", lw=0.5, ls=(0, (3, 2)), zorder=4)
    ax.text(n * 0.99, track.cut, f"{(1 - track.cut) * 100:.0f}%", ha="right", va="bottom",
            fontsize=5, color="#444444")
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(0, top * 1.06)
    ax.set_xticks([])
    ax.set_ylabel("distance", fontsize=5.5, color="#444444", labelpad=2)
    ax.tick_params(axis="y", labelsize=5, width=0.5, color="black", length=2)
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(4))
    for name, sp in ax.spines.items():
        sp.set_visible(name == "left")
        sp.set_linewidth(0.5)
    if title:
        ax.set_title(title, fontsize=5, linespacing=1.6)

    paths = []
    for fmt in formats:
        path = f"{out_stem}.dendrogram.{fmt}"
        fig.savefig(path, format=fmt)
        paths.append(path)
    plt.close(fig)
    return paths


def draw(track, out_stem, panel_mm=45.0, formats=("png", "pdf"), title=None,
         method="dendrogram"):
    """Draw the monomer tree. `method` is "dendrogram" (the grouping tree) or "nj"."""
    if method != "nj":
        return draw_dendrogram(track, out_stem, panel_mm, formats, title)
    return draw_nj(track, out_stem, panel_mm, formats, title)


def draw_nj(track, out_stem, panel_mm=45.0, formats=("png", "pdf"), title=None):
    """Write <out_stem>.tree.png/.pdf. Returns the paths, or [] if there is no tree."""
    ident = track.identity
    if ident is None or ident.shape[0] < 3:
        return []
    full = [u for u in track.units if not u.partial and u.satellite]
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
