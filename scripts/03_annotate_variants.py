#!/usr/bin/env python
"""Turn somatic variants into mutant proteins.

For every variant that falls in the coding part of a canonical transcript I
work out where it lands in the coding sequence, apply it, translate, and
compare the result with the normal protein. The output says what the variant
does to the protein and where the first changed residue is, which is what the
peptide step needs.

Things worth knowing about the implementation:

* Variants on the minus strand are applied in coding orientation, so the
  reference and alternative alleles are reverse-complemented first.
* The reference allele in the VCF has to match the reference base in the
  coding sequence. If it does not, the variant is counted and skipped rather
  than silently producing a wrong peptide.
* Indels that run off the end of a coding block (into an intron or a UTR) are
  skipped, because the protein consequence then depends on splicing.
* A variant that destroys the stop codon is reported but produces no peptides:
  the new C-terminus is read from the 3' UTR, which this pipeline does not load.

Usage:
  python scripts/03_annotate_variants.py <ref_dir> <out_dir> <vcf.gz> [<vcf.gz> ...]
"""
import gzip
import pickle
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import revcomp, translate, write_table

BUCKET = 100_000  # for finding the transcripts near a position


def build_position_index(transcripts):
    """Map (chromosome, 100 kb bucket) to the transcripts overlapping it."""
    index = defaultdict(list)
    for transcript, record in transcripts.items():
        chrom = record["chrom"]
        for start, end in record["blocks"]:
            for bucket in range(start // BUCKET, end // BUCKET + 1):
                index[(chrom, bucket)].append(transcript)
    return {key: sorted(set(value)) for key, value in index.items()}


def cds_index(record, position):
    """Position of a chromosomal coordinate in the coding sequence, or None."""
    offset = 0
    for start, end in record["blocks"]:
        if start <= position <= end:
            if record["strand"] == "+":
                return offset + (position - start)
            return offset + (end - position)
        offset += end - start + 1
    return None


def read_vcf(path):
    """Yield PASS variants as (chrom, pos, ref, alt, vaf)."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or "PASS" not in fields[6]:
                continue
            chrom, pos, ref, alt, info = fields[0], int(fields[1]), fields[3], fields[4], fields[7]
            if "," in alt or alt in (".", "") or not set(ref + alt) <= set("ACGTacgt"):
                continue  # multi-allelic or symbolic alleles are out of scope
            vaf = ""
            for entry in info.split(";"):
                if entry.startswith("TVAF="):
                    vaf = entry[5:]
                    break
            yield chrom.removeprefix("chr"), pos, ref.upper(), alt.upper(), vaf


def apply_variant(record, pos, ref, alt):
    """Return the mutant coding sequence, or a reason why the variant is skipped."""
    span = list(range(pos, pos + len(ref)))
    indices = [cds_index(record, p) for p in span]
    if any(i is None for i in indices):
        return None, "outside the coding blocks"

    if record["strand"] == "+":
        start = indices[0]
        cds_ref, cds_alt = ref, alt
    else:
        start = indices[-1]
        cds_ref, cds_alt = revcomp(ref), revcomp(alt)

    if sorted(indices) != list(range(start, start + len(ref))):
        return None, "spans a coding block boundary"

    cds = record["cds"]
    if cds[start:start + len(cds_ref)].upper() != cds_ref:
        return None, "reference allele does not match the reference sequence"

    return cds[:start] + cds_alt + cds[start + len(cds_ref):], None


def classify(wt_protein, mt_protein, length_change, hit_stop):
    """Name the consequence and find the first changed residue (1-based)."""
    if mt_protein == wt_protein:
        return "synonymous", None

    first = next(
        (i for i, (a, b) in enumerate(zip(wt_protein, mt_protein)) if a != b),
        min(len(wt_protein), len(mt_protein)),
    )

    if length_change % 3:
        return "frameshift", first + 1
    if not hit_stop and len(mt_protein) >= len(wt_protein):
        # The stop codon itself was changed, so translation ran to the end.
        return "stop_lost", first + 1
    if length_change == 0:
        if len(mt_protein) < len(wt_protein):
            return "stop_gained", first + 1
        changed = sum(1 for a, b in zip(wt_protein, mt_protein) if a != b)
        return ("missense" if changed == 1 else "multiple_substitutions"), first + 1
    if length_change > 0:
        return "inframe_insertion", first + 1
    return "inframe_deletion", first + 1


def main(ref_dir, out_dir, *vcfs):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    with gzip.open(Path(ref_dir) / "transcripts.pkl.gz", "rb") as handle:
        transcripts = pickle.load(handle)
    index = build_position_index(transcripts)
    print(f"Loaded {len(transcripts)} transcripts")

    counts = defaultdict(int)
    rows = []
    for vcf in vcfs:
        seen = 0
        for chrom, pos, ref, alt, vaf in read_vcf(vcf):
            seen += 1
            candidates = index.get((chrom, pos // BUCKET), [])
            if not candidates:
                counts["variants outside coding sequence"] += 1
                continue
            overlaps_coding = False
            for transcript in candidates:
                record = transcripts[transcript]
                if cds_index(record, pos) is None:
                    continue
                overlaps_coding = True
                mutant_cds, problem = apply_variant(record, pos, ref, alt)
                if problem:
                    counts[f"skipped: {problem}"] += 1
                    continue
                mt_protein, hit_stop = translate(mutant_cds)
                wt_protein = record["protein"]
                consequence, first = classify(
                    wt_protein, mt_protein, len(alt) - len(ref), hit_stop
                )
                counts[consequence] += 1
                change = ""
                if consequence in ("missense", "stop_gained") and first:
                    change = (
                        f"p.{wt_protein[first - 1]}{first}"
                        f"{mt_protein[first - 1] if first <= len(mt_protein) else '*'}"
                    )
                rows.append({
                    "variant_id": f"{chrom}:{pos}:{ref}>{alt}",
                    "chrom": chrom, "pos": pos, "ref": ref, "alt": alt, "vaf": vaf,
                    "transcript": transcript,
                    "gene_id": record["gene_id"], "gene_name": record["gene_name"],
                    "consequence": consequence,
                    "protein_change": change,
                    "first_changed_residue": first or "",
                    "wt_length": len(wt_protein), "mt_length": len(mt_protein),
                    "wt_protein": wt_protein, "mt_protein": mt_protein,
                })
            if not overlaps_coding:
                counts["variants not annotated on any transcript"] += 1
        print(f"  {Path(vcf).name}: {seen} PASS variants read")

    write_table(out / "variants.tsv", rows)

    with open(out / "mutant_proteins.faa", "w") as handle:
        for row in rows:
            if row["consequence"] == "synonymous":
                continue
            handle.write(
                f">{row['variant_id']}|{row['transcript']}|{row['gene_name']}"
                f"|{row['consequence']}\n{row['mt_protein']}\n"
            )

    with open(out / "annotation_summary.tsv", "w") as handle:
        handle.write("item\tcount\n")
        for key, value in sorted(counts.items(), key=lambda kv: -kv[1]):
            handle.write(f"{key}\t{value}\n")

    print(f"{len(rows)} variant/transcript pairs annotated")
    for key, value in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {key}: {value}")
    print(f"Done. Tables are in {out}/")


if __name__ == "__main__":
    main(*sys.argv[1:])
