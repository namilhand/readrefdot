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
With `--monomer`, also `<readid>.tree.png` / `.tree.pdf`.

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

**Colour.** **Black** (`--colour_main`) is the *main diagonal*: the alignment's own
diagonals plus the self-identity one, each then followed along its full continuous
extent through collinear matches. **Grey70** (`--colour_ext`) is every other diagonal —
the repeat ladder. Reverse-complement matches are vermillion (an inversion or foldback)
and are never a main diagonal.

Following the diagonal past the aligned part is the point rather than a side effect.
When an insertion duplicates the reference just upstream, the duplicated copy sits on
the same diagonal as the alignment, so black traces both copies and you can read
straight off the plot how much of the reference the read carries twice. The alignment
is used to find these diagonals; it is never drawn.

**Axes.** Ticks every 5 kb in both blocks, labelled only at each block's two ends —
at 45 mm there is no room for more, and the ends give you the range. Reference
coordinates are Mb to 3 decimals (`6.981`), read coordinates kb to 1 (`21.5`). The read
axis starts at 0 unless the primary alignment carries a **leading hard clip**, in which
case it starts at that offset (soft-clipped bases are present in SEQ and shift nothing).

The plot box is `--panel-mm` square (45 mm by default), excluding title and labels. The
reference window is exactly the read's aligned span, with no padding.

## How the main diagonal is determined

Black is not "the longest line" — it is derived from the alignment, in three steps.

**1. Seed a band per aligned block.** Each alignment on the chosen chromosome (primary
plus any supplementary in the same orientation) is walked through its CIGAR. Every
aligned block of at least `MIN_BLOCK` (20 bp) yields a *band*: the diagonal it sits on
plus the x-range it covers. Each block contributes two bands, one on the reference side
of the plot and its mirror on the read side. The self-identity diagonal is added as a
band spanning the whole plot.

**2. Grow each band along its own diagonal.** A band is extended through k-mer runs that
lie on the same diagonal (within `DIAG_TOL`, 3 bp) and are separated from it by at most
`EXTEND_GAP` (200 bp), repeating until it stops growing. This is what makes a black line
a continuous diagonal instead of only the aligned part of one, and it is what puts a
tandem duplication in black: the duplicated copy sits on the alignment's own diagonal,
just displaced.

**3. Colour each run.** A drawn run is black when its diagonal is within `DIAG_TOL` of a
band's and it overlaps that band's x-range by at least `MIN_OVERLAP` (half the run's own
length). Requiring real overlap rather than any overlap keeps a run that merely clips a
band's edge, or sits at the same offset kilobases away, in the ladder colour.

`MIN_BLOCK`, `DIAG_TOL`, `EXTEND_GAP` and `MIN_OVERLAP` are constants at the top of
`plot.py`, not options. `EXTEND_GAP` is the one worth knowing about: it sets how far a
black line may jump to keep following its diagonal. Raise it and black follows through
longer interruptions; in a very dense repeat array a large value could chain further
than you intend.

How far the black path reaches is otherwise governed by `--min-seg`: if an alignment
ends in a block shorter than `--min-seg`, no run there survives the length filter and
the line stops that far short of the panel edge.

## Options

| option | default | notes |
|---|---|---|
| `-k, --kmer` | 20 | seed size (5–31). A random 20-mer match has probability ~4⁻²⁰ |
| `--min-seg` | 170 | drop diagonal runs shorter than this (bp). **The main de-cluttering control** — 170 is just under one CEN178 satellite unit (units vary 177–179 bp), so each surviving run is at least one monomer and the ladder spacing reads as the repeat period |
| `--merge-gap` | k+1 | largest gap chained into one run. k+1 is exactly the step across a single substitution, so runs bridge isolated SNPs and nothing more; `1` gives strictly exact runs |
| `--panel-mm` | 45 | size of the square plot box in mm, excluding title and labels |
| `--outdir` | `.` | output directory |
| `--name` | read id | output file stem, giving `<name>.quad.png` / `.pdf`. One read only |
| `--colour_main` | black | colour of the main diagonals (alignment diagonals, followed to their full extent) |
| `--colour_ext` | grey70 | colour of every other diagonal |
| `--monomer` | off | annotate satellite monomers instead of coordinates (below) |
| `--monomer-period` | detect | satellite unit length in bp |
| `--monomer-cut` | 0.95 | identity at which two monomers join the same group |
| `--monomer-style` | lollipop | `lollipop` (stick + circle) or `block` |
| `--monomer-consensus` | built-in CEN178 | FASTA of the repeat consensus that fixes where a unit starts; `none` derives the phase from the sequence |
| `--monomer-tsv` | off | also write `<name>.monomers.tsv`, one row per unit |
| `--no-tree` | off | skip the monomer tree |

## Monomer annotation (`--monomer`)

For a centromeric read, genomic coordinates say little: both the read and the reference
window are tandem arrays of the ~178 bp CEN178 (aTha178) satellite unit. `--monomer`
replaces the coordinate axes with a strip of coloured blocks along the **top and the
right of the whole panel** — one block per monomer, coloured by similarity group — so
the axes say what the array is *made of*.

```bash
readrefdot --bam sample.bam --ref Col-0.fa --read "…/85266687/ccs" --monomer
```

Because the strips run the full concatenated axis, they annotate the reference block and
the read block in turn, and the two can be read against each other: an insertion of
whole satellite units shows up as extra blocks repeating the reference's colour pattern.

Nothing is taken from the aligner. With a consensus (the default) each unit is placed by
its own alignment to that consensus; without one, the units are found in the sequences
themselves. Four steps:

1. **Period.** Almost every k-mer in a tandem array recurs one unit later, so the
   histogram of distances between successive copies of the same 16-mer has a sharp mode
   at the unit length. No mode ⇒ not an array ⇒ the plot falls back to coordinates.
2. **Phase — from the published consensus by default.** Where a unit *starts* is the one
   thing the sequence cannot settle on its own: any rotation of the monomer tiles it
   equally well, and which rotation you get depends on which self-anchor happens to rank
   first, so two reads over the same locus can be cut at different points. So the phase is
   taken from the published CEN178 consensus, which ships with the tool. Every consensus
   k-mer that matches at sequence position `b` from consensus position `a` implies a unit
   start at `b − a`; those implied starts are the votes. Orientation is detected (arrays
   carry the monomer either way round — Chr4 reverse, Chr1 and Chr3 forward in Col-0), and
   if the consensus does not phase the array the tool falls back to the self-anchored
   method below. `--monomer-consensus FASTA` supplies a different repeat's consensus;
   `--monomer-consensus none` forces the self-anchored method.

   The two agree: on four reads, consensus tiling and self-anchored tiling produce the
   same units differing by a **constant rotation** (18, 6, 7 and 86 bp), each exactly
   `178 −` the rotation measured independently between the array consensus and the
   published one.

   **The votes only say where to start; each unit is then placed by its own alignment.**
   A semi-global Needleman–Wunsch puts the whole consensus inside a one-unit window, and
   the next window starts where that alignment ended. This is what vote offsets alone
   cannot do: an indel *inside* a unit splits that unit's k-mers into two clusters (those
   before it vote at one position, those after it at another) and a tolerance can only
   pick one of them. An alignment absorbs it — a 2 bp insertion 10 bp into a unit makes
   that unit 180 bp and leaves all 79 downstream units at 178 bp with unchanged identity.

   The window is deliberately one unit wide plus slack: a wider one lets the alignment
   pick whichever unit ahead scores best and silently skip the one in front of it (which
   marked every other unit non-satellite when first tried). The gap penalty must exceed
   the mismatch penalty for the same reason it does in any aligner — at 1/−1/−1 a fifth of
   the units came out 177 bp because a gap explained a substitution as cheaply as a
   mismatch; at 2/−3/−5 the same units are 178 bp at identical identity.

   **Non-satellite stretches are labelled, not tiled.** A window the consensus does not
   match at ≥60% identity is emitted as a non-satellite block rather than cut into
   pretend units: a synthetic 1,200 bp insert of random sequence comes back as two blocks
   totalling 1,025 bp, and real centromeric reads have none. Non-satellite blocks are dark
   grey on the axis, excluded from the groups and from the tree, and carry `satellite=0`
   in `--monomer-tsv`; every satellite unit carries its `identity` to the consensus
   (median 0.93–0.96 on real reads, and as low as 0.74 for genuinely degenerate units).

3. **Phase without a consensus (the fallback).** Every 16-mer that recurs at the period is
   a candidate marker of the same point in successive units. Each is checked for a *consistent* offset from the best one
   along the whole array — modulo the period, so a marker still votes in the units where
   the best one was mutated away — and a k-mer sitting at a different place in different
   units is not a phase marker and is dropped. The survivors (up to 60) vote on where a
   unit starts. **No single marker survives in every unit**: the best one covers 86–97% of
   them. The panel as a whole does — on real reads every unit contains at least one marker,
   median 43–53 of them, and 98–99% of boundaries land on a vote backed by dozens of
   agreeing markers. Reference and read are tiled separately but from the same anchor
   panel, so a boundary means the same thing in both blocks.

   **The votes place the boundaries; gaps are filled, not walked.** Every unit carries its
   own votes, so boundaries are read straight off them, and a gap of about *m* periods gets
   *m*−1 boundaries spread evenly across it (~1% of boundaries, never more than one in a
   row). An earlier version walked one period at a time and snapped to a vote only within a
   quarter period, which could not survive an indel: after an insertion of L bp every
   downstream vote sits L off the walk's targets, and when `L mod 178` fell outside the snap
   window (45–133 bp) the walk never met a vote again and placed the whole rest of the block
   blind and out of phase. Reading the votes directly keeps an indel local — the unit
   containing it comes out long, and the next unit starts where its own votes say.
4. **Groups.** Units are compared by shared 8-mer content (a Mash-style identity
   estimate — exact pairwise alignment would be more accurate and is not worth it, since
   the numbers only have to separate satellite variants) and clustered by average
   linkage, cut at `--monomer-cut`. Groups are colour-ordered by size, so the most
   abundant variant is always the first palette colour. Partial units at the array edges
   are left light grey, non-satellite blocks dark grey.

Each unit is drawn as a **lollipop** — a short stick from the panel edge with a circle
centred on its end, coloured by group. The strip is only as thick as those marks need,
so the annotation sits tight against the panel. The circle is auto-sized to the space
one unit actually gets: a 45 mm
panel over ~220 monomers leaves each 0.2 mm, so at that density the circles sit side by
side and the strip reads as a coloured line on a comb rather than as separate dots.
`--monomer-style block` draws a solid block per unit instead, which is denser to read at
a glance and is the better choice for a long array.

Grouping is **per plot**: a colour identifies a variant within one figure and carries no
meaning across figures. `--monomer-tsv` writes the units out — position, group, length,
sequence — which is what to use when groups need to be compared between plots.

The default cut of 0.95 is what separates CEN178 variants in these arrays; 0.90 merges
almost everything into one group, and above 0.97 the groups start to split on individual
substitutions. The title reports how many units, of what length, in how many groups.

### The tree

`--monomer` also writes `<name>.tree.png` / `.pdf`: a neighbour-joining tree of the same
units, tips coloured by the same groups, laid out with Felsenstein's equal-angle
algorithm so it reads as an unrooted radial tree. The dot plot says where the units are;
the tree says how they are related. `--no-tree` skips it.

Two things to know when reading it. The distance is **alignment-free** — 1 − the
shared-k-mer identity estimate that drives the grouping — so the scale bar is labelled
"distance", not substitutions per site; it is a similarity tree of satellite variants,
not a substitution-model phylogeny. And **identical units land on the same point**, so a
tree of 200 monomers can show 40 visible tips; that overlap is real information (those
monomers have the same sequence), not a drawing fault.

### Relation to the published CEN178 consensus

The units this finds are the published CEN178 monomer, but the tool derives its own
phase from the sequence, so boundaries do not start where the published consensus
starts. Rotated into register, a per-array consensus built from the tiled units matches
the published 178 bp consensus at **94–98%** identity:

| array | identity | rotation | strand |
|---|---|---|---|
| Chr4 Col-0 | 94.4% (10/178 differ) | 160 bp | reverse |
| Chr4 Ler-0 | 97.8% (4/178) | 172 bp | reverse |
| Chr1 Col-0 | 95.5% (8/178) | 171 bp | forward |

Since the tool now takes its phase from that consensus by default, unit 1 starts at the
consensus start in every plot, and unit coordinates in `--monomer-tsv` are directly
comparable between reads, samples and papers. The `phase` column records which reference
fixed them, and the orientation the array carries.

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
| `monomer`, `monomer-period`, `monomer-cut`, `monomer-style`, `monomer-consensus` | | monomer annotation; `monomer` is on for anything but `0`/`no`/`false`. Rows with it on also write `<suffix>.tree.png`/`.pdf` |

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

## contrib

`contrib/` holds helpers that assume a particular upstream data format and are not part
of the tool — currently `fill_annotation.py`, which fills a manifest's annotation
columns from CHARLA de-novo indel tables. See `contrib/README.md`.

## Requirements

Python ≥3.8, `pysam`, `numpy`, `matplotlib`.
