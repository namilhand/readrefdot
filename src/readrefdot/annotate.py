"""Optional guide lines at coordinates the user supplies."""

from dataclasses import dataclass, field
from typing import List


def parse_positions(text):
    """'6988767,6992861' or '6,988,767 6992861' -> [6988767, 6992861]."""
    if not text:
        return []
    out = []
    for tok in str(text).replace(";", ",").replace(" ", ",").split(","):
        tok = tok.replace(",", "").strip()
        if tok:
            out.append(int(tok))
    return out


def pair_up(points):
    """Endpoints -> intervals.

    The lists hold DISTINCT endpoints, so how they pair depends on the count. Four values
    are two separate intervals and pair off two by two. But an insertion sitting
    immediately beside its donor -- a tandem duplication -- shares an endpoint with it, so
    the two intervals arrive as three values and have to be read as a chain. Pairing those
    two by two would keep the insertion and silently drop the donor."""
    pts = sorted(set(points))
    z = zip(pts[::2], pts[1::2]) if len(pts) % 2 == 0 else zip(pts[:-1], pts[1:])
    return [(a, b) for a, b in z if b > a]


@dataclass
class Lines:
    """`ref` are 1-based absolute reference positions; `read` are 0-based read positions."""
    ref: List[int] = field(default_factory=list)
    read: List[int] = field(default_factory=list)

    @classmethod
    def parse(cls, ref_lines=None, read_lines=None):
        return cls(ref=parse_positions(ref_lines), read=parse_positions(read_lines))

    def __bool__(self):
        return bool(self.ref or self.read)

    def intervals(self, ctx):
        """(reference intervals, read intervals) on the concatenated [ref | read] axis.

        The same coordinates `marks` returns, paired into the segments they delimit and
        clipped to their own block; an interval entirely outside the plotted window is
        dropped, exactly as a guide line for it would be."""
        R, Q = len(ctx.ref_seq), len(ctx.read_seq)
        out = []
        for pts, base, lo, hi in ((self.ref, -1 - ctx.win_start, 0, R),
                                  (self.read, R, R, R + Q)):
            got = []
            for a, b in pair_up(p + base for p in pts):
                a, b = max(lo, a), min(hi, b)
                if b > a:
                    got.append((a, b))
            out.append(got)
        return out[0], out[1]

    def marks(self, ctx):
        """Positions on the concatenated [reference | read] axis of the quad plot.

        Each coordinate is clipped to its OWN block. A reference position outside the
        plotted window maps past the block boundary, and would otherwise be drawn inside
        the read block -- a line in the wrong half of the plot rather than no line."""
        R, Q = len(ctx.ref_seq), len(ctx.read_seq)
        out = []
        for p in self.ref:
            v = p - 1 - ctx.win_start
            if 0 <= v <= R:
                out.append(v)
        for p in self.read:
            if 0 <= p <= Q:
                out.append(R + p)
        return out
