# contrib

Helpers that ship with readrefdot but are **not part of the tool**. They assume a
particular upstream data format, so they are useful only if you produce that format.
readrefdot itself has no dependency on anything here.

## `fill_annotation.py` — CHARLA de-novo indels

Fills the `ref-lines` / `read-lines` columns of a `readrefdot-batch` manifest from
CHARLA's de-novo indel tables, so insertions, deletions and INS donor regions are marked
on the plots.

```bash
python3 contrib/fill_annotation.py manifest.tsv \
    --all-tsv   .../02-indels/output/all_denovo.tsv \
    --donor-tsv .../02-indels/output/donor_tracing_denovo.tsv
readrefdot-batch manifest.tsv --force
```

| event | `ref-lines` | `read-lines` |
|---|---|---|
| INS | donor region on the reference | inserted sequence + the donor's copy in the read |
| DEL | deleted block on the reference | deletion point in the read |

The script exists because three details are easy to get wrong by hand:

**Read coordinates are mirrored for minus-strand reads.** CHARLA stores them in FASTQ
orientation; readrefdot draws the read in alignment (SEQ) orientation, so
`seq = read_length - fastq`.

**The donor's position in the read is deduced, and the sign flips with strand.** It is
not tabulated. The offset is (gap between donor and INS site) + insertion size, applied
to the right of the insertion when the donor is right of the INS site *in the
reference* — a rule stated in reference orientation, so it inverts on a minus-strand
read. Every result is checked against a projection of the donor through the alignment
CIGAR, and flagged `CHECK:` if the two disagree by more than 20 bp.

**Identity comes from the BAM filename**, `<part_id>.<hap>.sorted.bam`, not from the
manifest's `suffix` — a suffix can name a different part than its own BAM.

Coordinates outside the plotted window, or past the end of the read, are still written:
the manifest stays a complete record and readrefdot clips each coordinate to its own
half of the plot at draw time. A distant (non-tandem) donor legitimately produces
`donor not covered by this alignment`.
