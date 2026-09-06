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

    def marks(self, ctx):
        """Positions on the concatenated [reference | read] axis of the quad plot."""
        R = len(ctx.ref_seq)
        n = R + len(ctx.read_seq)
        vals = [p - 1 - ctx.win_start for p in self.ref] + [R + p for p in self.read]
        return [v for v in vals if 0 <= v <= n]
