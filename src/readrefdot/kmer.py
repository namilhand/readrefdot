"""Exact k-mer matching engine: sequence -> diagonal segments.

Nothing here knows about BAMs, plots or files -- it turns two strings into runs of
collinear exact k-mer matches, which is all a dot plot is.
"""

import numpy as np

COMP = bytes.maketrans(b"ACGTacgtNn", b"TGCAtgcaNn")

_CODE = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate("ACGT"):
    _CODE[ord(_b)] = _i
    _CODE[ord(_b.lower())] = _i


def revcomp(s):
    return s.translate(COMP)[::-1]


def kmer_codes(seq, k):
    """Encode every k-mer as a uint64 (k <= 31). Returns (values, valid_mask)."""
    c = _CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    n = len(seq) - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.uint64), np.empty(0, dtype=bool)
    bad = (c == 255).astype(np.int32)
    cum = np.concatenate(([0], np.cumsum(bad)))
    valid = (cum[k:] - cum[:-k]) == 0
    c64 = c.astype(np.uint64)
    vals = np.zeros(n, dtype=np.uint64)
    for j in range(k):
        vals = (vals << np.uint64(2)) | c64[j:j + n]
    return vals, valid


def build_index(seq, k, mask=None):
    """Sorted k-mer index of `seq` for match(). `mask` may invalidate positions."""
    vals, valid = kmer_codes(seq, k)
    if mask is not None:
        valid = valid & mask
    order = np.flatnonzero(valid)
    sorted_vals = vals[order]
    o = np.argsort(sorted_vals, kind="stable")
    return vals, valid, sorted_vals[o], order[o]


def match(query_vals, query_valid, sorted_vals, order, max_hits=0):
    """All exact (query_pos, target_pos) k-mer matches, fully vectorised.

    max_hits <= 0 means no repeat masking. Returns (qpos, tpos, n_skipped)."""
    if max_hits <= 0:
        max_hits = np.iinfo(np.int64).max
    qpos = np.flatnonzero(query_valid)
    if qpos.size == 0:
        return np.empty(0, np.int64), np.empty(0, np.int64), 0
    qv = query_vals[qpos]
    lo = np.searchsorted(sorted_vals, qv, "left")
    hi = np.searchsorted(sorted_vals, qv, "right")
    cnt = hi - lo
    n_skipped = int(np.count_nonzero(cnt > max_hits))
    keep = (cnt > 0) & (cnt <= max_hits)
    lo, cnt, qpos = lo[keep], cnt[keep], qpos[keep]
    total = int(cnt.sum())
    if total == 0:
        return np.empty(0, np.int64), np.empty(0, np.int64), n_skipped
    off = np.cumsum(cnt) - cnt
    idx = np.repeat(lo, cnt) + (np.arange(total) - np.repeat(off, cnt))
    return np.repeat(qpos, cnt), order[idx].astype(np.int64), n_skipped


def merge_segments(a_pos, b_pos, gap, reverse=False):
    """Chain collinear k-mer hits into runs. Returns (a0, b0, a1, b1) arrays.

    Forward runs share `b - a`; reverse-complement runs share `b + a`. `gap` is the
    largest step in `a` that still continues a run -- it bridges single SNPs (which
    kill k consecutive k-mers) and holes left by repeat masking."""
    if a_pos.size == 0:
        z = np.empty(0, np.int64)
        return z, z, z, z
    key = (b_pos + a_pos) if reverse else (b_pos - a_pos)
    order = np.lexsort((a_pos, key))
    key, a, b = key[order], a_pos[order], b_pos[order]
    brk = np.empty(a.size, dtype=bool)
    brk[0] = True
    brk[1:] = (key[1:] != key[:-1]) | (np.diff(a) > gap)
    s = np.flatnonzero(brk)
    e = np.append(s[1:], a.size) - 1
    return a[s], b[s], a[e], b[e]


def filter_min_length(seg, k, min_seg):
    """Drop runs shorter than `min_seg` bp. The main de-cluttering control."""
    if min_seg <= 0 or seg[0].size == 0:
        return seg
    keep = (seg[2] - seg[0] + k) >= min_seg
    return tuple(arr[keep] for arr in seg)


def compare(query, target, k, merge_gap=None, min_seg=0, max_hits=0,
            target_mask=None, query_mask=None):
    """Full pipeline: two sequences -> (forward_runs, reverse_runs, stats).

    Runs are (a0, b0, a1, b1) with `a` indexing `query` and `b` indexing `target`."""
    if merge_gap is None:
        merge_gap = k + 1
    _, _, sorted_vals, order = build_index(target, k, target_mask)

    qv, qvalid = kmer_codes(query, k)
    if query_mask is not None:
        qvalid = qvalid & query_mask
    fa, fb, skip_f = match(qv, qvalid, sorted_vals, order, max_hits)

    rc = revcomp(query)
    cv, cvalid = kmer_codes(rc, k)
    if query_mask is not None:
        cvalid = cvalid & query_mask[::-1]
    cq, rb, skip_r = match(cv, cvalid, sorted_vals, order, max_hits)
    ra = (len(query) - k - cq).astype(np.int64)

    fwd = filter_min_length(merge_segments(fa, fb, merge_gap, False), k, min_seg)
    rev = filter_min_length(merge_segments(ra, rb, merge_gap, True), k, min_seg)
    stats = dict(n_fwd_hits=int(fa.size), n_rev_hits=int(ra.size),
                 n_skipped=int(skip_f + skip_r), merge_gap=merge_gap)
    return fwd, rev, stats
