#!/usr/bin/env python
"""Compare the pipeline's output on the simulated data with the planted truth.

Prints one line per check and exits with a non-zero status if any of them
fails, so the test can be used as a gate before trusting a real run.

Usage: python scripts/check_test_results.py <test_data_dir> <results_dir>
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import read_table

PASS, FAIL = "ok  ", "FAIL"


class Checks:
    def __init__(self):
        self.failed = 0

    def __call__(self, condition, description, detail=""):
        print(f"  [{PASS if condition else FAIL}] {description}"
              + (f" -- {detail}" if detail and not condition else ""))
        if not condition:
            self.failed += 1
        return condition


def main(test_dir, results_dir):
    truth = json.loads((Path(test_dir) / "truth.json").read_text())
    results = Path(results_dir)

    variants = read_table(results / "variants.tsv")
    peptides = read_table(results / "peptides.tsv")
    candidates = read_table(results / "candidates.tsv")
    annotation = {row["item"]: int(row["count"])
                  for row in read_table(results / "annotation_summary.tsv")}
    peptide_counts = {row["item"]: int(row["count"])
                      for row in read_table(results / "peptide_summary.tsv")}

    by_variant = {row["variant_id"]: row for row in variants}
    check = Checks()

    print("Annotation")
    for expected in truth["variants"]:
        row = by_variant.get(expected["variant_id"])
        if not check(row is not None,
                     f"{expected['gene']}: variant was annotated",
                     expected["variant_id"]):
            continue
        check(row["consequence"] == expected["consequence"],
              f"{expected['gene']}: called {expected['consequence']}",
              f"got {row['consequence']}")
        check(row["gene_name"] == expected["gene"],
              f"{expected['gene']}: assigned to the right gene",
              f"got {row['gene_name']}")
        if expected.get("protein_change"):
            check(row["protein_change"] == expected["protein_change"],
                  f"{expected['gene']}: protein change {expected['protein_change']}",
                  f"got {row['protein_change'] or 'nothing'}")

    print("Variants that have to be rejected")
    for entry in truth["rejected"]:
        check(entry["variant_id"] not in by_variant,
              f"not annotated: {entry['reason']}",
              f"{entry['variant_id']} was annotated anyway")
        key = next((k for k in annotation if entry["reason"] in k), None)
        check(key is not None and annotation[key] >= 1,
              f"reported as skipped: {entry['reason']}",
              "no such line in annotation_summary.tsv")

    print("Peptides")
    peptides_by_variant = {}
    for row in peptides:
        peptides_by_variant.setdefault(row["variant_id"], []).append(row)

    for variant_id in truth["no_peptides_from"]:
        check(variant_id not in peptides_by_variant,
              f"no peptides from {by_variant.get(variant_id, {}).get('consequence', variant_id)}",
              f"{len(peptides_by_variant.get(variant_id, []))} peptides")

    missense = [e for e in truth["variants"] if e["consequence"] == "missense"]
    for expected in missense:
        rows = peptides_by_variant.get(expected["variant_id"], [])
        check(len(rows) > 0, f"{expected['gene']}: produced peptides")
        lengths = {int(r["length"]) for r in rows}
        check(lengths <= {8, 9, 10, 11} and 9 in lengths,
              f"{expected['gene']}: peptide lengths are 8 to 11",
              f"got {sorted(lengths)}")
        differences = [
            sum(1 for a, b in zip(r["peptide"], r["wt_peptide"]) if a != b)
            for r in rows if r["wt_peptide"]
        ]
        check(differences and set(differences) == {1},
              f"{expected['gene']}: every peptide differs from the normal one "
              "at exactly one position",
              f"got {sorted(set(differences))}")

    frameshift = next(e for e in truth["variants"] if e["consequence"] == "frameshift")
    rows = peptides_by_variant.get(frameshift["variant_id"], [])
    check(len(rows) > 0, "frameshift: produced peptides")
    check(all(not r["wt_peptide"] for r in rows),
          "frameshift: no normal counterpart is claimed")

    deletion = next(e for e in truth["variants"]
                    if e["consequence"] == "inframe_deletion")
    check(len(peptides_by_variant.get(deletion["variant_id"], [])) > 0,
          "in-frame deletion: produced peptides")

    print("Self filter")
    decoy = truth["self_filtered_peptide"]
    check(decoy not in {r["peptide"] for r in peptides},
          f"the peptide that also occurs in a normal protein ({decoy}) was removed")
    check(peptide_counts.get("dropped: also present in a normal protein", 0) >= 1,
          "the self filter reported at least one removal")

    print("Expression filter")
    silent_gene = truth["gene_not_expressed"]
    check(any(r["gene_name"] == silent_gene for r in peptides),
          f"{silent_gene}: peptides were generated")
    check(not any(r["gene_name"] == silent_gene for r in candidates),
          f"{silent_gene}: dropped because the gene is not expressed")

    print("Binding predictions (published HLA-A*02:01 epitopes as positive controls)")
    ranks = {}
    for row in candidates:
        ranks.setdefault(row["peptide"], int(row["rank"]))
    for epitope in truth["positive_controls"]:
        if check(epitope in ranks, f"{epitope} is among the candidates"):
            rank = ranks[epitope]
            check(rank <= 5, f"{epitope} is ranked in the top 5", f"ranked {rank}")

    print()
    if check.failed:
        print(f"{check.failed} check(s) failed.")
        return 1
    print("All checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main(*sys.argv[1:3]))
