#!/usr/bin/env python
"""Build a small transcript database for annotating coding variants.

For every canonical protein-coding transcript I store where its coding bases
sit on the chromosome, together with its coding sequence. That is all that is
needed to turn a genomic variant into a mutant protein.

Two details matter here:

* Ensembl's GTF lists CDS features without the stop codon, while the CDS FASTA
  includes it. I therefore read both the `CDS` and the `stop_codon` features.
* A transcript is only kept if the coding blocks on the chromosome add up to
  the length of its coding sequence, and if translating that sequence gives
  back exactly the protein Ensembl distributes. Transcripts that fail (for
  example ones with a ribosomal slippage or an incomplete 3' end) are dropped,
  because a silent coordinate shift there would invent mutant peptides that do
  not exist.

Usage:
  python scripts/02_build_reference.py <gtf.gz> <cds.fa.gz> <pep.fa.gz> <out_dir>
"""
import gzip
import pickle
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import read_fasta, strip_version, translate

ATTRIBUTE = re.compile(r'(\S+) "([^"]*)"')


def parse_gtf(path):
    """Collect coding blocks per canonical protein-coding transcript."""
    blocks = defaultdict(list)
    info = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in ("CDS", "stop_codon"):
                continue
            attributes = dict(ATTRIBUTE.findall(fields[8]))
            if attributes.get("transcript_biotype") != "protein_coding":
                continue
            if "Ensembl_canonical" not in fields[8]:
                continue
            transcript = strip_version(attributes["transcript_id"])
            blocks[transcript].append((int(fields[3]), int(fields[4])))
            info[transcript] = {
                "chrom": fields[0],
                "strand": fields[6],
                "gene_id": strip_version(attributes.get("gene_id", "")),
                "gene_name": attributes.get("gene_name", ""),
            }
    return blocks, info


def main(gtf_path, cds_path, pep_path, out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)

    print("Reading the GTF")
    blocks, info = parse_gtf(gtf_path)
    print(f"  {len(blocks)} canonical protein-coding transcripts")

    print("Reading coding sequences")
    cds = {strip_version(header.split()[0]): seq
           for header, seq in read_fasta(cds_path)}

    # The protein FASTA is keyed by protein id, but its headers name the
    # transcript each protein comes from, which is how I match the two up.
    print("Reading proteins")
    pep_by_tx = {}
    opener = gzip.open if str(pep_path).endswith(".gz") else open
    with opener(pep_path, "rt") as handle:
        current = None
        for line in handle:
            if line.startswith(">"):
                match = re.search(r"transcript:(\S+)", line)
                current = strip_version(match.group(1)) if match else None
                if current:
                    pep_by_tx[current] = []
            elif current:
                pep_by_tx[current].append(line.strip())
    pep_by_tx = {k: "".join(v).rstrip("*") for k, v in pep_by_tx.items()}

    print("Checking each transcript against its own coding sequence")
    transcripts = {}
    dropped = defaultdict(int)

    for transcript, block_list in blocks.items():
        meta = info[transcript]
        sequence = cds.get(transcript)
        if not sequence:
            dropped["no coding sequence in the FASTA"] += 1
            continue
        ordered = sorted(block_list, reverse=meta["strand"] == "-")
        span = sum(end - start + 1 for start, end in ordered)
        if span != len(sequence):
            dropped["coding blocks and sequence have different lengths"] += 1
            continue
        if len(sequence) % 3:
            dropped["coding sequence is not a multiple of three"] += 1
            continue
        protein, _ = translate(sequence)
        expected = pep_by_tx.get(transcript)
        if expected is None:
            dropped["no protein in the FASTA"] += 1
            continue
        if protein != expected:
            dropped["translation does not match the Ensembl protein"] += 1
            continue
        transcripts[transcript] = {
            "chrom": meta["chrom"],
            "strand": meta["strand"],
            "gene_id": meta["gene_id"],
            "gene_name": meta["gene_name"],
            "blocks": ordered,
            "cds": sequence,
            "protein": protein,
        }

    print(f"  kept {len(transcripts)}")
    for reason, count in sorted(dropped.items(), key=lambda kv: -kv[1]):
        print(f"  dropped {count}: {reason}")

    # Every protein Ensembl knows about, isoforms included, as one string.
    # A mutant peptide that occurs anywhere in here is not tumour-specific.
    print("Writing the normal proteome for the self filter")
    with open(out / "self_proteome.txt", "w") as handle:
        handle.write("|".join(pep_by_tx.values()))

    with gzip.open(out / "transcripts.pkl.gz", "wb") as handle:
        pickle.dump(transcripts, handle, protocol=4)

    with open(out / "reference_summary.tsv", "w") as handle:
        handle.write("item\tvalue\n")
        handle.write(f"canonical transcripts in the GTF\t{len(blocks)}\n")
        handle.write(f"transcripts kept\t{len(transcripts)}\n")
        for reason, count in sorted(dropped.items()):
            handle.write(f"dropped: {reason}\t{count}\n")
    print(f"Done. Reference is in {out}/")


if __name__ == "__main__":
    main(*sys.argv[1:5])
