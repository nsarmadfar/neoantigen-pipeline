#!/usr/bin/env python
"""Predict how well each peptide is presented by each of the sample's HLA alleles.

MHCflurry 2 is used here. Its presentation score combines predicted binding
affinity with an antigen-processing model, and the presentation percentile puts
that score on a scale where low numbers mean "unusually good for this allele",
which is what makes scores comparable between alleles.

Every peptide is scored against every allele separately, rather than only
against the best one, because which allele presents a peptide matters for
interpretation, and because an allele-level summary is useful on its own.

The matching normal peptides are scored too. A mutant peptide that binds much
better than its normal counterpart is the more promising candidate, since the
immune system has not been trained to tolerate it.

Usage:
  python scripts/05_predict_binding.py <hla.tsv> <results_dir>
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import read_table


def load_alleles(path):
    alleles = []
    for line in Path(path).read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        first = line.split("\t")[0].strip()
        if first.lower() == "allele":
            continue
        alleles.append(first)
    return alleles


def load_predictor():
    try:
        from mhcflurry import Class1PresentationPredictor
    except ImportError:
        sys.exit("MHCflurry is not installed. conda env create -f environment.yml")
    try:
        return Class1PresentationPredictor.load()
    except Exception as error:  # the models are a separate download
        sys.exit(
            f"Could not load the MHCflurry models ({error}).\n"
            "Run: mhcflurry-downloads fetch models_class1_presentation"
        )


def predict(predictor, peptides, alleles):
    """Score every peptide against every allele, one allele at a time."""
    frames = []
    for allele in alleles:
        print(f"  {allele}: {len(peptides)} peptides")
        result = predictor.predict(
            peptides=peptides,
            alleles={"sample": [allele]},
            include_affinity_percentile=True,
            verbose=0,
        )
        result = result[[
            "peptide", "affinity", "affinity_percentile",
            "processing_score", "presentation_score", "presentation_percentile",
        ]].copy()
        result.insert(1, "allele", allele)
        frames.append(result)
    return pd.concat(frames, ignore_index=True)


def main(hla_path, results_dir):
    results = Path(results_dir)
    alleles = load_alleles(hla_path)
    print(f"HLA alleles: {', '.join(alleles)}")

    predictor = load_predictor()
    supported = set(predictor.supported_alleles)
    unsupported = [a for a in alleles if a not in supported]
    if unsupported:
        sys.exit(
            f"MHCflurry has no model for: {', '.join(unsupported)}.\n"
            "Remove the allele from config/hla.tsv or use another predictor for it."
        )

    rows = read_table(results / "peptides.tsv")
    if not rows:
        print("No peptides to score.")
        for name in ("binding_mutant.tsv", "binding_normal.tsv"):
            (results / name).write_text(
                "peptide\tallele\taffinity\taffinity_percentile\t"
                "processing_score\tpresentation_score\tpresentation_percentile\n"
            )
        return

    mutant = sorted({row["peptide"] for row in rows})
    normal = sorted({row["wt_peptide"] for row in rows if row["wt_peptide"]})

    print(f"Scoring {len(mutant)} mutant peptides")
    predict(predictor, mutant, alleles).to_csv(
        results / "binding_mutant.tsv", sep="\t", index=False
    )

    print(f"Scoring {len(normal)} matching normal peptides")
    predict(predictor, normal, alleles).to_csv(
        results / "binding_normal.tsv", sep="\t", index=False
    )

    print(f"Done. Predictions are in {results}/")


if __name__ == "__main__":
    main(*sys.argv[1:3])
