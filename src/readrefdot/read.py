"""BAM access: a read's sequence and the reference window it maps to.

The BAM is used for exactly two things -- the read's bases, and the coordinates that
say which slice of the reference to compare against. Nothing the aligner decided about
*how* the read aligns is carried into the plot.
"""

import os
from dataclasses import dataclass, field
from typing import List

import pysam

from .kmer import revcomp

FLANK = 0             # bp of reference padding around the aligned span


class ReadNotFound(Exception):
    pass


@dataclass
class ReadContext:
    read_id: str
    chrom: str
    win_start: int                 # 0-based inclusive
    win_end: int                   # 0-based exclusive
    ref_seq: str
    read_seq: str                  # alignment orientation (revcomp if the read is reverse)
    strand: str
    aln_start: int = 0             # 0-based reference start of the read's alignment
    aln_end: int = 0               # 0-based exclusive reference end
    read_offset: int = 0           # read coord of read_seq[0]: leading hard clip, if any
    full_read_len: int = 0         # read length including hard-clipped bases
    alignments: List[object] = field(default_factory=list)   # on `chrom`
    primary: object = None
    n_alignments: int = 1
    n_secondary_ignored: int = 0

    @property
    def window(self):
        return f"{self.chrom}:{self.win_start + 1:,}-{self.win_end:,}"

    @property
    def read_len(self):
        return len(self.read_seq)


def open_fasta(path):
    """pysam.FastaFile, tolerating a symlink whose .fai sits beside the real file."""
    for p in (path, os.path.realpath(path)):
        if os.path.exists(p + ".fai"):
            return pysam.FastaFile(p)
    return pysam.FastaFile(path)


def _is_usable(a):
    return not (a.is_unmapped or a.is_secondary)


def iter_primary(bam_path):
    """Yield every primary alignment in the BAM, in file order."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for a in bam.fetch(until_eof=True):
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue
            if a.query_sequence:
                yield a.query_name


def _context_from(alns, n_sec, fasta, read_id):
    primary = next((a for a in alns if not a.is_supplementary), alns[0])
    if primary.query_sequence is None:
        raise ReadNotFound(f"{read_id}: primary alignment carries no SEQ")

    span = {}
    for a in alns:
        span[a.reference_name] = span.get(a.reference_name, 0) + a.reference_length
    chrom = max(span, key=span.get)
    use = [a for a in alns if a.reference_name == chrom]

    aln_s = min(a.reference_start for a in use)
    aln_e = max(a.reference_end for a in use)
    win_s = max(0, aln_s - FLANK)
    win_e = min(fasta.get_reference_length(chrom), aln_e + FLANK)

    seq = primary.query_sequence          # BAM SEQ is already alignment orientation
    # SEQ omits hard-clipped bases, so read_seq[0] is not read position 0 when the
    # primary carries a leading hard clip. Soft clips are present in SEQ and add nothing.
    ct = primary.cigartuples or []
    offset = ct[0][1] if ct and ct[0][0] == 5 else 0
    full_len = sum(l for op, l in ct if op in (0, 1, 4, 5, 7, 8))
    return ReadContext(read_id=read_id, chrom=chrom, win_start=win_s, win_end=win_e,
                       ref_seq=fasta.fetch(chrom, win_s, win_e), read_seq=seq,
                       strand="-" if primary.is_reverse else "+",
                       aln_start=aln_s, aln_end=aln_e, read_offset=offset, full_read_len=full_len or len(seq),
                       alignments=use, primary=primary,
                       n_alignments=len(alns), n_secondary_ignored=n_sec)


def load(bam_path, ref_path, read_ids, skip_missing=False):
    """Yield a ReadContext for each requested read id, in one pass over the BAM.

    With skip_missing=True a read that is absent yields a ReadNotFound instance (carrying
    .read_id) instead of raising, so one missing read cannot abort the rest of the batch.

    Secondary alignments (flag 0x100) are alternative placements, so they are ignored:
    including them would stretch the reference window for no gain."""
    wanted = set(read_ids)
    found, n_sec = {}, {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for a in bam.fetch(until_eof=True):
            if a.query_name not in wanted or a.is_unmapped:
                continue
            if a.is_secondary:
                n_sec[a.query_name] = n_sec.get(a.query_name, 0) + 1
                continue
            found.setdefault(a.query_name, []).append(a)

    fasta = open_fasta(ref_path)
    try:
        for rid in read_ids:
            if rid not in found:
                err = ReadNotFound(f"read {rid!r} not found in {bam_path}")
                err.read_id = rid
                if skip_missing:
                    yield err
                    continue
                raise err
            yield _context_from(found[rid], n_sec.get(rid, 0), fasta, rid)
    finally:
        fasta.close()
