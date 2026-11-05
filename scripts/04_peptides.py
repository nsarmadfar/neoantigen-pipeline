#!/usr/bin/env python
"""Cut the mutant proteins into the peptides that MHC class I could present.

Class I peptides are mostly 8 to 11 residues long, so for every protein change
I take all windows of those lengths that cover the changed part of the protein.
For a substitution that is the mutated residue; for a frameshift it is
everything from the shift to the new stop codon.

Every peptide is then checked against the whole normal proteome, isoforms
included. A mutant peptide that also exists somewhere in a normal protein is
not tumour-specific, and the immune system is tolerant to it, so it is dropped.
This also removes windows that happen to fall next to the mutation but do not
actually contain the change.

Stop-gained and stop-lost variants are reported in the variant table but
produce no peptides here: a premature stop only shortens the protein, and the
new C-terminus of a stop-lost protein is read from the 3' UTR, which this
pipeline does not load.

Usage:
  python scripts/04_peptides.py <ref_dir> <results_dir> <lengths, e.g. 8,9,10,11>
"""
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import AMINO_ACIDS, read_table, write_table

# Consequences that can create a peptide the normal proteome does not contain.
PEPTIDE_YIELDING = {
    "missense", "multiple_substitutions",
    "inframe_insertion", "inframe_deletion", "frameshift",
}


def changed_region(row):
    """Residue positions (0-based, in the mutant protein) that are new."""
    wt, mt = row["wt_protein"], row["mt_protein"]
    first = int(row["first_changed_residue"]) - 1
    if row["consequence"] == "frameshift":
        return first, len(mt) - 1
    if row["consequence"] in ("missense", "multiple_substitutions"):
        last = max(
            (i for i, (a, b) in enumerate(zip(wt, mt)) if a != b),
            default=first,
        )
        return first, last
    # In-frame indel: the inserted residues, or the junction left by a deletion.
    inserted = max(0, len(mt) - len(wt))
    return first, first + inserted


def windows(protein, first, last, lengths):
    """All windows of the given lengths that overlap residues first..last."""
    for length in lengths:
        lowest = max(0, first - length + 1)
        highest = min(last, len(protein) - length)
        for start in range(lowest, highest + 1):
            yield start, protein[start:start + length]


def find_in_normal_proteome(proteome, wanted):
    """Which of our peptides also occur in a normal protein.

    The obvious way round is to ask, for each peptide, whether it is in the
    proteome. That means scanning 48 MB once per peptide. Walking the proteome
    once and looking each window up in a set of our own peptides does the same
    job in one pass, which is thousands of times less work.
    """
    found = set()
    for length, peptides in sorted(wanted.items()):
        if not peptides:
            continue
        print(f"  scanning the proteome for {len(peptides)} peptides of length {length}")
        for start in range(len(proteome) - length + 1):
            window = proteome[start:start + length]
            if window in peptides:
                found.add(window)
    return found


def main(ref_dir, results_dir, lengths_arg):
    results = Path(results_dir)
    lengths = [int(x) for x in lengths_arg.split(",")]
    self_proteome = (Path(ref_dir) / "self_proteome.txt").read_text()
    print(f"Normal proteome for the self filter: {len(self_proteome):,} residues")

    rows = read_table(results / "variants.tsv")
    counts = defaultdict(int)
    candidates = []
    wanted = defaultdict(set)
    seen = set()

    for row in rows:
        if row["consequence"] not in PEPTIDE_YIELDING:
            continue
        if not row["first_changed_residue"]:
            continue
        first, last = changed_region(row)
        mt = row["mt_protein"]
        wt = row["wt_protein"]
        same_length = len(mt) == len(wt)
        for start, peptide in windows(mt, first, last, lengths):
            counts["windows over the change"] += 1
            if not set(peptide) <= AMINO_ACIDS:
                counts["dropped: unusual residue (X, U or a stop)"] += 1
                continue
            key = (peptide, row["variant_id"], row["transcript"])
            if key in seen:
                continue
            seen.add(key)
            wanted[len(peptide)].add(peptide)
            candidates.append({
                "peptide": peptide,
                "length": len(peptide),
                # The same window of the normal protein, for comparing how well
                # the mutant and the normal version bind. Only defined when the
                # mutation did not shift the sequence.
                "wt_peptide": wt[start:start + len(peptide)] if same_length else "",
                "position_in_protein": start + 1,
                "variant_id": row["variant_id"],
                "transcript": row["transcript"],
                "gene_id": row["gene_id"],
                "gene_name": row["gene_name"],
                "consequence": row["consequence"],
                "protein_change": row["protein_change"],
                "vaf": row["vaf"],
            })

    print(f"{len(candidates)} windows over a protein change; now the self filter")
    in_normal = find_in_normal_proteome(self_proteome, wanted)
    peptides = [row for row in candidates if row["peptide"] not in in_normal]
    counts["dropped: also present in a normal protein"] = len(candidates) - len(peptides)
    counts["kept"] = len(peptides)

    write_table(results / "peptides.tsv", peptides)
    with open(results / "peptide_summary.tsv", "w") as handle:
        handle.write("item\tcount\n")
        for key, value in sorted(counts.items(), key=lambda kv: -kv[1]):
            handle.write(f"{key}\t{value}\n")

    print(f"{len(peptides)} mutant peptides from {len({p['variant_id'] for p in peptides})} variants")
    for key, value in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {key}: {value}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
