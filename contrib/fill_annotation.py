#!/usr/bin/env python3
"""Fill ref-lines / read-lines in a readrefdot manifest from CHARLA indel tables.

    python3 fill_annotation.py input/manifest.tsv

For each manifest row it looks the event up in CHARLA's de-novo tables and writes the
two annotation columns:

  INS   ref-lines  = donor region on the reference   (rep_start, rep_end)
        read-lines = inserted sequence + the donor's copy in the read
  DEL   ref-lines  = deleted block on the reference  (start, end)
        read-lines = deletion point in the read      (read_start, read_end)

Three things this handles that are easy to get wrong:

* Read coordinates. CHARLA stores them in FASTQ orientation; readrefdot draws the read
  in alignment (SEQ) orientation. For a minus-strand read the two are mirrored, so every
  read coordinate is converted with seq = read_length - fastq.

* The donor's position in the read. It is not tabulated, so it is deduced from the
  reference geometry: offset = (distance between the donor and the INS site) + size,
  applied to the right of the insertion when the donor is to the right in the reference
  and to the left when it is to the left. That sign is expressed in reference
  orientation, so it INVERTS for a minus-strand read. Each result is checked against a
  projection of the donor through the alignment CIGAR.

* Which event a row refers to. part_id and hap_region come from the BAM filename
  (<part_id>.<hap>.sorted.bam), never from the suffix -- a suffix can disagree with its
  own BAM, and the BAM is authoritative.

Coordinates that fall outside the plotted window or past the end of the read are still
written; readrefdot clips each to its own half of the plot and simply draws fewer lines.
"""

import argparse
import csv
import os
import re
import sys

import pysam

DONOR_TSV = ("/Users/namilhand/1_LocalWork/01_Cambridge/01_CHARLA/02_hifi_analysis/"
             "02-indels/output/donor_tracing_denovo.tsv")
ALL_TSV = ("/Users/namilhand/1_LocalWork/01_Cambridge/01_CHARLA/02_hifi_analysis/"
           "02-indels/output/all_denovo.tsv")
BAM_RE = re.compile(r"(?P<part>.+)\.(?P<hap>cc|lc|ca|la)\.sorted\.bam$")
EVENT_RE = re.compile(r"((?:INS|DEL|CONV)_\d+)")
PROJ_TOL = 20          # bp: how far the deduced donor may sit from the CIGAR projection


def load_tables(all_tsv, donor_tsv):
    events, by_read = {}, {}
    with open(all_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            events[(r["part_id"], r["hap_region"], r["id"])] = r
            by_read.setdefault(r["readid"], []).append(r)
    donors = {}
    if donor_tsv and os.path.exists(donor_tsv):
        with open(donor_tsv) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                donors[(r["part_id"], r["hap_region"], r["id"])] = r
    return events, by_read, donors


def find_event(row, events, by_read):
    """(part, hap) from the BAM filename; event id from the suffix. Falls back to the
    read id when the suffix and the BAM disagree."""
    m = BAM_RE.match(os.path.basename(row["bam"]))
    ev = EVENT_RE.search(row.get("suffix", "") or "")
    ev = ev.group(1) if ev else None
    if m and ev:
        e = events.get((m.group("part"), m.group("hap"), ev))
        if e is not None and e["readid"] == row["readid"]:
            return e, None
    cand = [r for r in by_read.get(row["readid"], []) if r["id"] == ev] or \
           by_read.get(row["readid"], [])
    if len(cand) == 1:
        note = "matched by read id (suffix/BAM disagree)" if m and ev else None
        return cand[0], note
    if not cand:
        return None, "no event found for this read id"
    return None, f"{len(cand)} events on this read; suffix did not disambiguate"


def primary_alignment(bam, read_id, chrom, pos):
    with pysam.AlignmentFile(bam, "rb") as b:
        for a in b.fetch(chrom, max(0, pos - 5), pos + 5):
            if (a.query_name == read_id and not a.is_supplementary
                    and not a.is_secondary and a.query_sequence):
                return a
    return None


def ref_to_read(aln, ref_pos):
    """Reference position -> index into the read's SEQ."""
    q, r = 0, aln.reference_start
    for op, ln in aln.cigartuples:
        if op in (0, 7, 8):
            if r <= ref_pos < r + ln:
                return q + (ref_pos - r)
            q += ln; r += ln
        elif op == 1:
            q += ln
        elif op in (2, 3):
            if r <= ref_pos < r + ln:
                return q
            r += ln
        elif op in (4, 5):
            q += ln
    return None


def lines_for(row, event, donor):
    """Return (ref_lines, read_lines, note)."""
    aln = primary_alignment(row["bam"], row["readid"], event["Chr"], int(event["start"]))
    if aln is None:
        return None, None, "read not found in its BAM"
    L = len(aln.query_sequence)
    rev = event["strand"] == "-"

    def to_seq(x, y):
        return (L - y, L - x) if rev else (x, y)

    rs, re_ = int(event["read_start"]), int(event["read_end"])
    ins = to_seq(rs, re_)

    if event["type"] != "INS" or not donor or donor.get("rep_chr") != event["Chr"]:
        why = "DEL" if event["type"] != "INS" else "no donor on this chromosome"
        return ([int(event["start"]), int(event["end"])], sorted(set(ins)), why)

    st, en = int(event["start"]), int(event["end"])
    ds, de = int(donor["rep_start"]), int(donor["rep_end"])
    size = int(event["size"])
    if ds > st:                                   # donor to the right in the reference
        side, dist = 1, max(0, ds - en)
    else:                                         # to the left, or overlapping
        side, dist = -1, max(0, st - de)
    off = (side if not rev else -side) * (dist + size)
    dread = to_seq(rs + off, re_ + off)

    proj = [v for v in (ref_to_read(aln, ds - 1), ref_to_read(aln, de - 1)) if v is not None]
    if len(proj) < 2:
        note = "donor not covered by this alignment"
    elif abs(dread[0] - min(proj)) > PROJ_TOL:
        note = f"CHECK: deduced donor {dread} vs CIGAR projection {tuple(sorted(proj))}"
    else:
        note = "ok"
    return [ds, de], sorted(set(list(ins) + list(dread))), note


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("manifest")
    p.add_argument("--donor-tsv", default=DONOR_TSV, help="CHARLA donor_tracing_*.tsv")
    p.add_argument("--all-tsv", default=ALL_TSV, help="CHARLA all_denovo.tsv")
    p.add_argument("-o", "--out", help="write here instead of updating in place")
    p.add_argument("--dry-run", action="store_true", help="report, write nothing")
    a = p.parse_args()

    for f in (a.manifest, a.all_tsv):
        if not os.path.exists(f):
            sys.exit(f"not found: {f}")
    events, by_read, donors = load_tables(a.all_tsv, a.donor_tsv)

    with open(a.manifest, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit(f"{a.manifest}: no data rows")
    cols = list(rows[0].keys())
    for c in ("ref-lines", "read-lines"):
        if c not in cols:
            cols.append(c)

    n_ok = n_warn = n_fail = 0
    print(f"{'suffix':<24}{'ref-lines':<26}{'read-lines':<28}note")
    for r in rows:
        ev, why = find_event(r, events, by_read)
        if ev is None:
            print(f"{r.get('suffix','?'):<24}{'-':<26}{'-':<28}SKIPPED: {why}")
            n_fail += 1
            continue
        donor = donors.get((ev["part_id"], ev["hap_region"], ev["id"]))
        ref, read, note = lines_for(r, ev, donor)
        if ref is None:
            print(f"{r.get('suffix','?'):<24}{'-':<26}{'-':<28}SKIPPED: {note}")
            n_fail += 1
            continue
        r["ref-lines"] = ",".join(map(str, ref))
        r["read-lines"] = ",".join(map(str, read))
        if why:
            note = f"{note}; {why}"
        if note.startswith("CHECK"):
            n_warn += 1
        else:
            n_ok += 1
        print(f"{r.get('suffix','?'):<24}{r['ref-lines']:<26}{r['read-lines']:<28}{note}")

    if a.dry_run:
        print("\ndry run: nothing written")
        return
    out = a.out or a.manifest
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {out}  ({n_ok} filled"
          + (f", {n_warn} to check" if n_warn else "")
          + (f", {n_fail} skipped" if n_fail else "") + ")")
    if n_fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
