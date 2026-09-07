"""readrefdot - k-mer dot plots of a long read against the reference it maps to.

    from readrefdot import load, Params, quad
    for ctx in load("sample.bam", "genome.fa", ["m84227_.../85266687/ccs"]):
        quad(ctx, Params(), "out")
"""

__version__ = "0.1.0"

from .annotate import Lines                              # noqa: F401
from .monomer import MonomerTrack, Unit                  # noqa: F401
from .monomer import annotate as annotate_monomers       # noqa: F401
from .plot import Params, quad                           # noqa: F401
from .read import ReadContext, ReadNotFound, load        # noqa: F401

__all__ = ["load", "ReadContext", "ReadNotFound", "Lines", "Params", "quad",
           "annotate_monomers", "MonomerTrack", "Unit", "__version__"]
