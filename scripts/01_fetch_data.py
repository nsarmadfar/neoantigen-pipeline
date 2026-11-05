#!/usr/bin/env python
"""Download the reference, the somatic variants and the expression table.

Files that are already there are left alone, so the script can be re-run
cheaply. The expression matrix from Expression Atlas covers 675 cell lines and
is about half a gigabyte; only the HCC1395 column is kept, in a small file.

Usage: python scripts/01_fetch_data.py config/settings.yaml data/raw
"""
import csv
import gzip
import sys
import urllib.request
from pathlib import Path

import yaml


def download(url, dest):
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  have {dest.name}")
        return dest
    print(f"  getting {dest.name}")
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(url, timeout=120) as response, open(tmp, "wb") as out:
        while chunk := response.read(1 << 20):
            out.write(chunk)
    tmp.rename(dest)
    return dest


def mean_of_list(cell):
    """Expression Atlas writes one value per assay, separated by commas."""
    values = [v for v in cell.strip().split(",") if v.strip()]
    if not values:
        return 0.0
    return sum(float(v) for v in values) / len(values)


def extract_expression(matrix, column, dest):
    """Pull one cell line out of the Expression Atlas TPM matrix."""
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  have {dest.name}")
        return dest
    print(f"  taking column {column} out of {matrix.name}")
    with open(matrix) as handle, open(dest, "w") as out:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        index = header.index(column)
        out.write("gene_id\tgene_name\ttpm\n")
        for row in reader:
            tpm = mean_of_list(row[index]) if len(row) > index else 0.0
            out.write(f"{row[0]}\t{row[1]}\t{tpm:g}\n")
    return dest


def check_assembly(vcf, assembly):
    """The variant coordinates have to match the reference we annotate against."""
    with gzip.open(vcf, "rt") as handle:
        for line in handle:
            if not line.startswith("##"):
                break
            if line.startswith("##reference") and assembly not in line:
                sys.exit(
                    f"{vcf.name} is not on {assembly}: {line.strip()}\n"
                    "Annotating against the wrong assembly would shift every position."
                )
    print(f"  {vcf.name}: coordinates are on {assembly}")


def main(config_path, out_dir):
    config = yaml.safe_load(Path(config_path).read_text())
    urls = config["downloads"]
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("Reference (Ensembl)")
    for key in ("gtf", "cds", "pep"):
        download(urls[key], out / Path(urls[key]).name)

    print("Somatic variants (SEQC2 truth set)")
    for key in ("snv", "indel"):
        vcf = download(urls[key], out / Path(urls[key]).name)
        check_assembly(vcf, config["assembly"])

    print("Expression (Expression Atlas)")
    matrix = download(urls["expression"], out / Path(urls["expression"]).name)
    extract_expression(matrix, urls["expression_column"], out / "expression.tsv")

    print(f"Done. Input files are in {out}/")


if __name__ == "__main__":
    main(*sys.argv[1:3])
