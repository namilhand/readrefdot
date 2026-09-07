"""CEN178 monomer annotation: find the satellite units, group them by sequence.

A centromeric read and the reference window it maps to are both tandem arrays of the
~178 bp CEN178 (aTha178) unit. This module cuts both sequences into those units and
sorts the units into similarity groups, so the dot plot can be labelled with what the
array is made of rather than with genomic coordinates.

Everything is derived from the sequences themselves -- no consensus, no library, no
aligner. Three steps:

  1. period      the spacing that recurs between identical k-mers (~178)
  2. phase       k-mers that recur at that spacing vote on where a unit starts, so a
                 unit boundary is a real motif and not an arbitrary offset
  3. groups      units are compared by shared k-mer content and clustered

Reference and read are tiled separately but with the SAME anchor panel, so a unit
boundary means the same thing in both blocks and the two annotations are comparable.
"""

from dataclasses import dataclass, field
from typing import List

import numpy as np

from .kmer import build_index, kmer_codes, match, revcomp

ANCHOR_K = 16              # k for period detection and phase anchors: long enough to be
                           # near-unique within one unit, short enough to survive a SNP
SIM_K = 8                  # k for unit-vs-unit similarity (a unit only has ~170 k-mers)
PERIOD_LO, PERIOD_HI = 150, 210     # search window for the unit length
PERIOD_TOL = 6             # bp a gap may differ from the period and still count
VOTE_TOL = 6               # bp within which boundary votes are the same boundary
MIN_ANCHOR_HITS = 3        # occurrences before a k-mer can be an anchor
MAX_ANCHORS = 60
CONS_K = 12                # k for matching the consensus: survives ~5% divergence
MIN_CONS_COVER = 0.5       # consensus phasing needs a vote for this many expected units
MIN_FULL = 0.7             # a unit shorter than this x period is a partial (array edge)

# The published CEN178 (aTha178) consensus. Used only to fix the PHASE -- where a unit
# starts -- so that unit 1 means the same stretch of the monomer in every plot and in the
# literature. Arrays match it at 94-98% and some carry it in reverse, both of which the
# matching handles; nothing here assumes the array equals the consensus.
CEN178 = ("AGTATAAGAACTTAAACCGCAACCCGATCTTAAAAGCCTAAGTAGTGTTTCCTTGTTAGAAGACACAAAGCCAAAGACTCA"
          "TATGGACTTTGGCTACACCATGAAAGCTTTGAGAAGCAAGAAGAAGGTTGGTTAGTGTTTTGGAGTCGAATATGACTTGAT"
          "GTCATGTGTATGATTG")
DEFAULT_CUT = 0.95         # group units whose estimated identity is >= this

# Group colours, most abundant group first. Okabe-Ito first (colour-blind safe), then
# distinguishable extras; groups past the end of the palette fall back to UNSET.
PALETTE = ["#D55E00", "#0072B2", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
           "#7B3294", "#B2DF23", "#8B4513", "#00CED1", "#F0E442", "#FF69B4"]
UNSET = "#DDDDDD"          # a group past the palette, or a partial unit


@dataclass
class Unit:
    """One monomer. `start`/`end` index the sequence it was tiled from."""
    block: str             # "ref" or "read"
    start: int
    end: int
    group: int = -1        # -1 = ungrouped (partial units)
    partial: bool = False

    @property
    def length(self):
        return self.end - self.start


@dataclass
class MonomerTrack:
    period: int
    phase: str = "de novo"             # what fixed the unit boundaries
    units: List[Unit] = field(default_factory=list)
    n_groups: int = 0
    identity: np.ndarray = None        # pairwise identity of the grouped units

    def colours(self):
        return [UNSET if u.group < 0 or u.group >= len(PALETTE) else PALETTE[u.group]
                for u in self.units]

    @property
    def n_full(self):
        return sum(1 for u in self.units if not u.partial)


def _kmer_positions(seq, k):
    """(values, positions) of every valid k-mer, sorted by value then position."""
    vals, valid = kmer_codes(seq, k)
    pos = np.flatnonzero(valid)
    v = vals[pos]
    o = np.lexsort((pos, v))
    return v[o], pos[o]


def detect_period(seq, k=ANCHOR_K, lo=PERIOD_LO, hi=PERIOD_HI):
    """The dominant spacing between successive copies of the same k-mer.

    In a tandem array almost every k-mer recurs one unit later, so the histogram of
    those gaps has a sharp mode at the unit length. Returns 0 if the sequence has no
    such periodicity (i.e. it is not a satellite array)."""
    v, p = _kmer_positions(seq, k)
    if v.size < 2:
        return 0
    same = v[1:] == v[:-1]
    d = np.diff(p)[same]
    d = d[(d >= lo) & (d <= hi)]
    if d.size < 10:
        return 0
    hist = np.bincount(d, minlength=hi + 1).astype(float)
    smooth = hist.copy()                       # +-1 bp, so 177/178/179 vote together
    smooth[1:-1] += hist[:-2] + hist[2:]
    peak = int(np.argmax(smooth[lo:hi + 1]) + lo)
    # The smoothed peak can land one bp off the true mode (177/178/179 vote together and
    # ties go to the lower index), so snap back to the raw mode beside it.
    period = max(range(max(lo, peak - 1), min(hi, peak + 1) + 1), key=lambda c: hist[c])
    # A real array: the modal gap must dominate. Otherwise call it non-periodic.
    if hist[period] < 0.02 * d.size:
        return 0
    return period


def _anchor_panel(seq, period, k=ANCHOR_K, max_anchors=MAX_ANCHORS):
    """k-mers that recur one period apart, best first. Each is a list of positions.

    These are the motifs that mark the same point in successive units. Using many of
    them, rather than one, means a unit whose best anchor was mutated away is still
    given a boundary by another anchor."""
    v, p = _kmer_positions(seq, k)
    if v.size == 0:
        return []
    starts = np.flatnonzero(np.r_[True, v[1:] != v[:-1]])
    ends = np.r_[starts[1:], v.size]
    expect = len(seq) / period
    out = []
    for s, e in zip(starts, ends):
        if e - s < MIN_ANCHOR_HITS or e - s > 3 * expect:
            continue                            # too rare, or several copies per unit
        q = p[s:e]
        good = int(np.count_nonzero(np.abs(np.diff(q) - period) <= PERIOD_TOL))
        if good >= MIN_ANCHOR_HITS - 1:
            out.append((good, q))
    out.sort(key=lambda t: -t[0])
    return [q for _, q in out[:max_anchors]]


def _boundary_votes(anchors, period):
    """Every anchor, shifted onto the phase of the best one, votes for a boundary.

    Returns (positions, weights) with weights = how many anchors support that boundary.
    An anchor is only used if its offset from the reference anchor is consistent along
    the whole array; a k-mer that occurs at a different place in different units is not
    a phase marker and is dropped."""
    if not anchors:
        return np.empty(0, int), np.empty(0, int)
    ref = anchors[0]
    shifted = [ref]
    for q in anchors[1:]:
        if q.size < MIN_ANCHOR_HITS:
            continue
        i = np.clip(np.searchsorted(ref, q), 1, ref.size - 1)
        near = np.where(np.abs(q - ref[i - 1]) <= np.abs(q - ref[i]), ref[i - 1], ref[i])
        # Offsets are taken MODULO the period. In the units where the reference anchor
        # itself is missing, the nearest reference copy is a whole period away -- and
        # those are exactly the units that need another anchor's vote. Folding the offset
        # keeps them instead of discarding them for being "too far".
        ang = 2 * np.pi * (q - near) / period
        centre = np.arctan2(np.sin(ang).mean(), np.cos(ang).mean()) * period / (2 * np.pi)
        resid = (q - near - centre + period / 2) % period - period / 2
        if np.median(np.abs(resid)) > VOTE_TOL:
            continue                            # inconsistent phase -> not an anchor
        shifted.append(np.rint(q - centre - np.median(resid)).astype(np.int64))

    allp = np.sort(np.concatenate(shifted))
    allp = allp[allp >= 0]
    # collapse votes that are the same boundary
    brk = np.flatnonzero(np.r_[True, np.diff(allp) > VOTE_TOL])
    seg_end = np.r_[brk[1:], allp.size]
    pos = np.array([int(np.median(allp[s:e])) for s, e in zip(brk, seg_end)])
    return pos, seg_end - brk


def _consensus_votes(seq, cons, k=CONS_K):
    """Unit starts implied by consensus k-mer matches, best orientation first.

    A consensus k-mer at consensus position `a` matching the sequence at `b` says a unit
    starts at `b - a`. Every surviving k-mer of every unit votes, so the phase comes from
    the consensus rather than from whichever self-anchor happened to rank first -- which
    is the one thing the self-anchored tiling cannot keep stable between reads.

    Returns (votes, weights, strand, n_hits)."""
    _, _, sorted_vals, order = build_index(seq, k)
    best = (np.empty(0, int), np.empty(0, int), "+", 0)
    for strand, c in (("+", cons), ("-", revcomp(cons))):
        cv, cvalid = kmer_codes(c, k)
        qpos, tpos, _ = match(cv, cvalid, sorted_vals, order)
        if qpos.size == 0:
            continue
        starts = np.sort(tpos - qpos)
        brk = np.flatnonzero(np.r_[True, np.diff(starts) > VOTE_TOL])
        end = np.r_[brk[1:], starts.size]
        pos = np.array([int(np.median(starts[a:b])) for a, b in zip(brk, end)])
        wt = end - brk
        if int(wt.sum()) > best[3]:
            best = (pos, wt, strand, int(wt.sum()))
    return best


def _tile_block(votes, weights, b0, b1, period):
    """Walk the votes into a continuous tiling of [b0, b1).

    From the strongest vote in the block, step one period at a time in both directions,
    snapping to a nearby vote when there is one and interpolating when there is not. A
    missing anchor therefore costs phase accuracy in that one unit, not the tiling."""
    sel = (votes >= b0) & (votes < b1)
    v, w = votes[sel], weights[sel]
    if v.size == 0:
        edges = list(range(b0, b1, period))
    else:
        seed = int(v[np.argmax(w)])
        win = max(4, period // 4)
        edges = [seed]
        for step in (period, -period):
            cur = seed
            while True:
                target = cur + step
                if not (b0 - period < target < b1 + period):
                    break
                near = np.flatnonzero(np.abs(v - target) <= win)
                if near.size:
                    cur = int(v[near[np.argmax(w[near])]])
                else:
                    cur = target
                if cur <= b0 - 1 or cur >= b1:
                    break
                edges.append(cur)
        edges.sort()
    edges = [e for e in edges if b0 <= e <= b1]
    bounds = sorted(set([b0] + edges + [b1]))
    return [(s, e) for s, e in zip(bounds[:-1], bounds[1:]) if e > s]


def _unit_similarity(seqs, k=SIM_K):
    """Estimated pairwise identity from shared k-mer content (a Mash-style estimate).

    Exact alignment of every pair would be more accurate and is not worth it here: the
    numbers only have to separate satellite variants, and Jaccard over 8-mers does that
    on 178 bp units while staying instant."""
    n = len(seqs)
    if n == 0:
        return np.zeros((0, 0))
    index, rows = {}, []
    for s in seqs:
        vals, valid = kmer_codes(s, k)
        keys = set(int(x) for x in vals[valid])
        rows.append(keys)
        for key in keys:
            index.setdefault(key, len(index))
    M = np.zeros((n, len(index)), dtype=np.float32)
    for i, keys in enumerate(rows):
        M[i, [index[key] for key in keys]] = 1.0
    inter = M @ M.T
    size = np.diag(inter).copy()
    union = size[:, None] + size[None, :] - inter
    jac = np.divide(inter, union, out=np.zeros_like(inter), where=union > 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        ident = 1.0 + np.log(np.clip(2 * jac / (1 + jac), 1e-9, 1.0)) / k
    ident = np.clip(ident, 0.0, 1.0)
    np.fill_diagonal(ident, 1.0)
    return ident


def _average_linkage(dist, cut):
    """Average-linkage agglomerative clustering, stopping at `cut`.

    n is a few hundred at most, so the naive O(n^3) loop is instant and avoids a scipy
    dependency."""
    n = dist.shape[0]
    if n == 0:
        return np.empty(0, dtype=int)
    d = dist.astype(float).copy()
    np.fill_diagonal(d, np.inf)
    members = [[i] for i in range(n)]
    alive = np.ones(n, dtype=bool)
    sizes = np.ones(n)
    while alive.sum() > 1:
        sub = np.where(alive[:, None] & alive[None, :], d, np.inf)
        i, j = np.unravel_index(np.argmin(sub), sub.shape)
        if sub[i, j] > cut:
            break
        new = (sizes[i] * d[i] + sizes[j] * d[j]) / (sizes[i] + sizes[j])
        d[i] = new
        d[:, i] = new
        d[i, i] = np.inf
        sizes[i] += sizes[j]
        members[i] += members[j]
        alive[j] = False
        members[j] = []
    labels = np.full(n, -1, dtype=int)
    groups = sorted((m for m in members if m), key=len, reverse=True)
    for g, mem in enumerate(groups):
        labels[mem] = g
    return labels


def annotate(ctx, period=None, cut=DEFAULT_CUT, anchor_k=ANCHOR_K, sim_k=SIM_K,
             consensus=CEN178):
    """Tile the reference window and the read into monomers and group them.

    `consensus` fixes the phase (default: the published CEN178 monomer); pass None to
    take the phase from the sequence itself. Either way the boundaries are placed by the
    same vote-and-walk, so the tiling is equally robust -- what the consensus buys is that
    unit 1 starts at the same point in the monomer in every plot.

    Returns a MonomerTrack whose unit coordinates index the concatenated
    [reference | read] axis of the quad plot. Returns None when the sequences are not
    a tandem array."""
    R = len(ctx.ref_seq)
    S = ctx.ref_seq + ctx.read_seq
    blocks = (("ref", 0, R), ("read", R, len(S)))

    if not period:
        period = detect_period(S, anchor_k)
    if not period:
        period = len(consensus) if consensus else 0
    if not period:
        return None

    phase = "de novo"
    votes = weights = None
    if consensus:
        v, w, strand, _ = _consensus_votes(S, consensus)
        # Only trust it if it actually phases the array: enough units carry a vote.
        strong = int(np.count_nonzero(w >= 3))
        if strong >= MIN_CONS_COVER * len(S) / period:
            votes, weights = v, w
            phase = f"consensus ({'forward' if strand == '+' else 'reverse'})"
    if votes is None:
        anchors = _anchor_panel(S, period, anchor_k)
        votes, weights = _boundary_votes(anchors, period)

    units = []
    for name, b0, b1 in blocks:
        for s, e in _tile_block(votes, weights, b0, b1, period):
            units.append(Unit(block=name, start=s, end=e,
                              partial=(e - s) < MIN_FULL * period))

    full = [u for u in units if not u.partial]
    if full:            # report the length units actually have, not the seed period
        lens = np.array([u.length for u in full])
        period = int(np.bincount(lens).argmax())
    ident = _unit_similarity([S[u.start:u.end] for u in full], sim_k)
    labels = _average_linkage(1.0 - ident, 1.0 - cut)
    for u, g in zip(full, labels):
        u.group = int(g)
    return MonomerTrack(period=period, phase=phase, units=units,
                        n_groups=int(labels.max() + 1) if labels.size else 0,
                        identity=ident)


def write_tsv(track, ctx, path):
    """One row per monomer: where it is, which group it is in, its sequence."""
    S = ctx.ref_seq + ctx.read_seq
    R = len(ctx.ref_seq)
    with open(path, "w") as fh:
        fh.write("block\tindex\tplot_start\tplot_end\tblock_start\tref_pos\t"
                 "length\tgroup\tpartial\tphase\tsequence\n")
        n = {"ref": 0, "read": 0}
        for u in track.units:
            n[u.block] += 1
            bstart = u.start if u.block == "ref" else u.start - R
            refpos = (ctx.win_start + u.start + 1) if u.block == "ref" else "."
            fh.write(f"{u.block}\t{n[u.block]}\t{u.start}\t{u.end}\t{bstart}\t{refpos}\t"
                     f"{u.length}\t{u.group if u.group >= 0 else '.'}\t"
                     f"{int(u.partial)}\t{track.phase}\t{S[u.start:u.end]}\n")
    return path
