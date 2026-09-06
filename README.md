# readrefdot

k-mer dot plots comparing a **long read** with the **reference genome region it maps to**.

Everything drawn is an **exact k-mer match found in the two sequences**. The aligner's
alignment is never drawn — it is used only to pick the reference window and to colour
which diagonals it placed the read on.

```bash
pip install -e .

readrefdot --bam sample.bam --ref genome.fa --read "m84227_.../85266687/ccs"
```

## Input

| | |
|---|---|
| `--bam` | BAM containing the read. Must be indexed for `--all` to be quick |
| `--ref` | the reference FASTA the BAM was aligned to |
| `--read ID [ID ...]` | read id(s) to plot |
| `--all` | plot every primary alignment in the BAM instead |

`--read` and `--all` are mutually exclusive; one is required.

The read is always drawn in **alignment orientation** — a reverse-strand read is
reverse-complemented, so its main diagonal runs bottom-left to top-right whichever
strand it came from.

Secondary alignments (flag 0x100) are ignored: they are alternative placements, and
including them would stretch the reference window for no gain.

## Output

`<outdir>/<readid>.quad.png` and `.quad.pdf`, at 300 dpi with fonts embedded.
Characters that cannot appear in a filename (`/` in particular) become `_`.

## The plot

`[reference | read]` on **both** axes, from a single self-comparison of the two
concatenated sequences. That one comparison fills four quadrants:

| quadrant | axes | shows |
|---|---|---|
| bottom-left | reference × reference | the reference window's own repeat structure |
| top-right | read × read | the read's own repeat structure |
| top-left | reference (x) × read (y) | read vs reference |
| bottom-right | read (x) × reference (y) | the transpose |

Reading the two corner quadrants against each other is the point: if the read's repeat
ladder continues the reference's unbroken, the read preserved the repeat phase.

**Colour.** Forward matches are split by whether the aligner placed the read on that
diagonal: **black** (`--colour_main`) for the diagonals of its aligned blocks plus the
self-identity diagonal, **grey70** (`--colour_ext`) for every other diagonal — the
repeat ladder. Reverse-complement matches are vermillion (an inversion or foldback);
they are never a main diagonal. This is the only use made of the alignment.

**Axes.** Ticks every 5 kb in both blocks, labelled only at each block's two ends —
at 45 mm there is no room for more, and the ends give you the range. Reference
coordinates are Mb to 3 decimals (`6.981`), read coordinates kb to 1 (`21.5`). The read
axis starts at 0 unless the primary alignment carries a **leading hard clip**, in which
case it starts at that offset (soft-clipped bases are present in SEQ and shift nothing).

The plot box is 45 mm square, excluding title and labels; the reference window is the
read's aligned span padded by 1 kb either side. Both are constants, not options.

## Options

| option | default | notes |
|---|---|---|
| `-k, --kmer` | 20 | seed size (5–31). A random 20-mer match has probability ~4⁻²⁰ |
| `--min-seg` | 170 | drop diagonal runs shorter than this (bp). **The main de-cluttering control** — 170 is just under one CEN178 satellite unit (units vary 177–179 bp), so each surviving run is at least one monomer and the ladder spacing reads as the repeat period |
| `--merge-gap` | k+1 | largest gap chained into one run. k+1 is exactly the step across a single substitution, so runs bridge isolated SNPs and nothing more; `1` gives strictly exact runs |
| `--outdir` | `.` | output directory |
| `--colour_main` | black | colour of the diagonals the aligner placed the read on |
| `--colour_ext` | grey70 | colour of every other diagonal |

## Batch — `readrefdot-batch`

Run many plots from a TSV, one row per plot. Each row carries its own inputs, output
location and arguments, so one file can mix samples, references and settings.

```bash
readrefdot-batch manifest.tsv            # --dry-run to preview, --force to redraw
```

| column | required | meaning |
|---|---|---|
| `bam` | ✓ | BAM containing the read |
| `readid` | ✓ | read to plot |
| `reference` | ✓ | reference FASTA |
| `outdir` | ✓ | directory to write into (created if absent) |
| `suffix` | ✓ | output file stem → `<outdir>/<suffix>.quad.png` and `.pdf` |
| `k`, `min-seg`, `merge-gap`, `panel-mm` | | per-row overrides; blank = default |
| `colour_main`, `colour_ext` | | per-row colours |
| `ref-lines`, `read-lines` | | annotation positions, comma-separated |

Column names accept either `-` or `_`. Blank cells mean "use the default".

**Rows whose output already exists are skipped**, so the wrapper is safe to re-run after
adding rows or after a failure; `--force` redraws.

The manifest is validated before anything is drawn — every unusable row is reported at
once, rather than failing partway through. Once running, a row that fails (a read absent
from its BAM, say) is reported and the rest continue.

Rows sharing a BAM and reference are grouped so the BAM is scanned once per group rather
than once per row; the scan dominates runtime.

## Annotation (optional)

```bash
readrefdot --bam sample.bam --ref genome.fa --read "…/85266687/ccs" \
    --ref-lines 6988767,6992861 --read-lines 6796,10891,14985
```

Draws dotted blue guide lines beneath the data. `--ref-lines` are absolute reference
positions (1-based), `--read-lines` are read positions (0-based, in the orientation
drawn). Both describe one read, so they cannot be combined with `--all`.

## What a line means

Each k-mer hit is an exact match. A drawn line is a run of hits on the **same
diagonal**, so it never absorbs an indel — a 1 bp indel shifts the diagonal and starts a
new line. Within a line, gaps up to `--merge-gap` are bridged, so a line is *exact
except for isolated single substitutions*. Use `--merge-gap 1` if you want lines that
are perfect matches end to end.

## Python API

```python
from readrefdot import load, Params, quad

for ctx in load("sample.bam", "genome.fa", ["m84227_.../85266687/ccs"]):
    quad(ctx, Params(min_seg=170), f"out/{ctx.read_id.replace('/', '_')}")
```

`readrefdot.kmer.compare(query, target, k)` is the engine on its own: two strings in,
diagonal runs out — no BAM, no plot.

## Requirements

Python ≥3.8, `pysam`, `numpy`, `matplotlib`.
