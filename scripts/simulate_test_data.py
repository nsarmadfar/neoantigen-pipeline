#!/usr/bin/env python
"""Make a small test dataset where I already know what the answer should be.

The files written here have exactly the same formats as the real inputs (an
Ensembl-style GTF, coding and protein FASTA files, somatic VCFs and a table of
expression values), so the pipeline runs on them unchanged. That is the point:
the test exercises the real code, not a simplified copy of it.

Ten toy genes on a toy chromosome carry one planted variant each:

  GENE1   plus strand, expressed        missense
  GENE2   minus strand, expressed       missense (checks the strand handling)
  GENE3   plus strand, expressed        synonymous: must give no peptides
  GENE4   plus strand, NOT expressed    missense: must be dropped on expression
  GENE5   plus strand, expressed        1 bp insertion -> frameshift
  GENE6   minus strand, expressed       3 bp deletion -> in-frame deletion
  GENE7   plus strand, expressed        nonsense: must give no peptides
  GENE8   plus strand, expressed        missense whose peptide also occurs in
                                        GENE1's normal protein, so the self
                                        filter has to remove it
  GENE9   plus strand, expressed        missense creating NLVPMVATV
  GENE10  plus strand, expressed        missense creating GILGFVFTL

The last two are published HLA-A*02:01 epitopes (the CMV pp65 and influenza A
M1 peptides). They are not neoantigens; they are here as positive controls,
because they are known to be presented by that allele. If the binding step is
wired up correctly they should come out at the top of the ranking.

Three more variants check that bad input is rejected rather than annotated:
one with a reference allele that does not match the reference sequence, one in
an intron, and one deletion that runs from an exon into an intron.

Usage: python scripts/simulate_test_data.py <out_dir>
"""
import gzip
import json
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import AMINO_ACIDS, CODON_TABLE, COMPLEMENT, translate

CHROM = "T1"
SEED = 20261105
CODONS_FOR = {}
for codon, residue in CODON_TABLE.items():
    CODONS_FOR.setdefault(residue, []).append(codon)

GENES = [
    # name, strand, exon count, TPM, length in residues
    ("GENE1", "+", 3, 120.0, 150),
    ("GENE2", "-", 2, 80.0, 130),
    ("GENE3", "+", 2, 60.0, 110),
    ("GENE4", "+", 2, 0.0, 110),
    ("GENE5", "+", 3, 45.0, 140),
    ("GENE6", "-", 2, 30.0, 120),
    ("GENE7", "+", 2, 25.0, 110),
    ("GENE8", "+", 2, 55.0, 120),
    ("GENE9", "+", 2, 95.0, 120),
    ("GENE10", "+", 2, 70.0, 120),
]


def single_base_changes(codon):
    """Every codon reachable from this one by changing one base."""
    for position in range(3):
        for base in "ACGT":
            if base != codon[position]:
                yield position, base, codon[:position] + base + codon[position + 1:]


def codon_for_change(residue, target=None, synonymous=False):
    """Find a codon of `residue` that a single base change turns into `target`.

    Returns (codon, position changed, new base, resulting residue).
    """
    for codon in CODONS_FOR[residue]:
        for position, base, new_codon in single_base_changes(codon):
            new_residue = CODON_TABLE[new_codon]
            if synonymous and new_residue == residue:
                return codon, position, base, new_residue
            if target and new_residue == target:
                return codon, position, base, new_residue
            if target is None and not synonymous and new_residue not in (residue, "*"):
                return codon, position, base, new_residue
    raise ValueError(f"no single-base route from {residue} to {target}")


def residue_that_mutates_into(target):
    """A residue, and the codon change, that becomes `target` in one base."""
    for residue in sorted(AMINO_ACIDS):
        if residue == target:
            continue
        try:
            codon, position, base, _ = codon_for_change(residue, target=target)
        except ValueError:
            continue
        return residue, codon, position, base
    raise ValueError(f"nothing mutates into {target}")


class Gene:
    """A toy gene: a protein, its codons, and where they sit on the chromosome."""

    def __init__(self, name, strand, exons, tpm, protein, codons, start):
        self.name = name
        self.transcript = "ENSTTEST" + name.removeprefix("GENE").zfill(4)
        self.protein_id = self.transcript.replace("ENST", "ENSP")
        self.strand = strand
        self.tpm = tpm
        self.protein = protein
        self.cds = "".join(codons) + "TAA"
        self.blocks = self._layout(exons, start)
        self.positions = self._positions()

    def _layout(self, exons, start):
        """Split the coding sequence over exons, with introns in between."""
        total = len(self.cds)
        # Keep the last exon comfortably longer than the stop codon.
        sizes = []
        remaining = total
        for index in range(exons - 1):
            size = (total // exons) // 3 * 3
            sizes.append(size)
            remaining -= size
        sizes.append(remaining)
        blocks, cursor = [], start
        for size in sizes:
            blocks.append((cursor, cursor + size - 1))
            cursor += size + 400  # intron
        # Blocks are in chromosome order; on the minus strand the transcript
        # reads them the other way round.
        return blocks

    def _positions(self):
        """Chromosome position of every coding base, in transcript order."""
        ordered = sorted(self.blocks, reverse=self.strand == "-")
        positions = []
        for start, end in ordered:
            span = range(start, end + 1)
            positions.extend(span if self.strand == "+" else reversed(span))
        assert len(positions) == len(self.cds)
        return positions

    def genomic_base(self, index):
        """The base on the forward strand at coding position `index`."""
        base = self.cds[index]
        return base if self.strand == "+" else base.translate(COMPLEMENT)

    def gtf_lines(self):
        """CDS and stop_codon features, the way Ensembl writes them."""
        attributes = (
            f'gene_id "{self.name}.1"; transcript_id "{self.transcript}"; '
            f'gene_name "{self.name}"; gene_biotype "protein_coding"; '
            f'transcript_biotype "protein_coding"; tag "Ensembl_canonical";'
        )
        stop_positions = set(self.positions[-3:])
        lines = []
        for start, end in self.blocks:
            coding = [p for p in range(start, end + 1) if p not in stop_positions]
            if coding:
                lines.append((CHROM, "test", "CDS", min(coding), max(coding),
                              ".", self.strand, "0", attributes))
            in_block = sorted(p for p in stop_positions if start <= p <= end)
            if in_block:
                lines.append((CHROM, "test", "stop_codon", min(in_block),
                              max(in_block), ".", self.strand, "0", attributes))
        return ["\t".join(str(f) for f in line) + "\n" for line in lines]


def build_genes(rng):
    genes, planted, cursor = {}, {}, 10_000

    # Plant the sequence features the variants need, then fill the rest randomly.
    def random_protein(length):
        return "M" + "".join(rng.choice(sorted(AMINO_ACIDS)) for _ in range(length - 1))

    specs = {}
    for name, strand, exons, tpm, length in GENES:
        specs[name] = dict(strand=strand, exons=exons, tpm=tpm,
                           protein=random_protein(length))

    # GENE8's peptide has to exist in GENE1's normal protein.
    gene1 = specs["GENE1"]["protein"]
    shared = gene1[60:69]
    decoy_middle = shared[4]
    stand_in, _, _, _ = residue_that_mutates_into(decoy_middle)
    specs["GENE8"]["protein"] = (
        specs["GENE8"]["protein"][:40]
        + shared[:4] + stand_in + shared[5:]
        + specs["GENE8"]["protein"][49:]
    )

    # The two positive controls, planted one residue away from the real epitope.
    for name, epitope in (("GENE9", "NLVPMVATV"), ("GENE10", "GILGFVFTL")):
        target = epitope[-1]
        stand_in, _, _, _ = residue_that_mutates_into(target)
        protein = specs[name]["protein"]
        specs[name]["protein"] = (
            protein[:50] + epitope[:-1] + stand_in + protein[59:]
        )

    # GENE7 needs a residue that one base change turns into a stop codon.
    specs["GENE7"]["protein"] = (
        specs["GENE7"]["protein"][:60] + "Q" + specs["GENE7"]["protein"][61:]
    )

    # Back-translate, choosing the codon each planted variant needs.
    changes = {
        "GENE1": ("missense", 70, None),
        "GENE2": ("missense", 55, None),
        "GENE3": ("synonymous", 40, None),
        "GENE4": ("missense", 45, None),
        "GENE5": ("frameshift", 60, None),
        "GENE6": ("inframe_deletion", 50, None),
        "GENE7": ("nonsense", 60, "*"),
        "GENE8": ("missense", 44, decoy_middle),
        "GENE9": ("missense", 58, "V"),
        "GENE10": ("missense", 58, "L"),
    }

    for name, strand, exons, tpm, _ in GENES:
        spec = specs[name]
        protein = spec["protein"]
        kind, residue_index, target = changes[name]
        codons = []
        for index, residue in enumerate(protein):
            if index == residue_index and kind in ("missense", "synonymous", "nonsense"):
                codon, position, base, result = codon_for_change(
                    residue,
                    target=target if kind in ("missense", "nonsense") else None,
                    synonymous=kind == "synonymous",
                )
                codons.append(codon)
                planted[name] = dict(kind=kind, residue_index=index,
                                     codon_offset=position, new_base=base,
                                     from_residue=residue, to_residue=result)
            else:
                codons.append(rng.choice(CODONS_FOR[residue]))
        if kind in ("frameshift", "inframe_deletion"):
            planted[name] = dict(kind=kind, residue_index=residue_index)

        gene = Gene(name, strand, exons, tpm, protein, codons, cursor)
        # Sanity: the coding sequence must translate back to the protein.
        assert translate(gene.cds)[0] == protein, name
        genes[name] = gene
        cursor = gene.blocks[-1][1] + 5_000

    return genes, planted


def make_variants(genes, planted):
    """Turn the planted changes into VCF records and the expected results."""
    snvs, indels, truth = [], [], []

    for name, change in planted.items():
        gene = genes[name]
        kind = change["kind"]

        if kind in ("missense", "synonymous", "nonsense"):
            index = change["residue_index"] * 3 + change["codon_offset"]
            pos = gene.positions[index]
            ref = gene.genomic_base(index)
            alt = (change["new_base"] if gene.strand == "+"
                   else change["new_base"].translate(COMPLEMENT))
            snvs.append((pos, ref, alt))
            consequence = {"missense": "missense", "synonymous": "synonymous",
                           "nonsense": "stop_gained"}[kind]
            entry = dict(
                variant_id=f"{CHROM}:{pos}:{ref}>{alt}", gene=name,
                consequence=consequence,
            )
            # A synonymous variant leaves the protein alone, so there is no
            # protein change to report for it.
            if consequence in ("missense", "stop_gained"):
                entry["protein_change"] = (
                    f"p.{change['from_residue']}{change['residue_index'] + 1}"
                    f"{change['to_residue']}")
            truth.append(entry)

        elif kind == "frameshift":
            index = change["residue_index"] * 3
            pos = gene.positions[index]
            ref = gene.genomic_base(index)
            indels.append((pos, ref, ref + "T"))
            truth.append(dict(variant_id=f"{CHROM}:{pos}:{ref}>{ref}T",
                              gene=name, consequence="frameshift"))

        elif kind == "inframe_deletion":
            # Delete one codon. On the minus strand the chromosome positions of
            # that codon run backwards, and the VCF anchor base is the one just
            # before the deletion in chromosome order.
            index = change["residue_index"] * 3
            deleted = [index, index + 1, index + 2]
            positions = sorted(gene.positions[i] for i in deleted)
            anchor = positions[0] - 1
            anchor_index = gene.positions.index(anchor)
            ref = "".join(
                gene.genomic_base(gene.positions.index(p))
                for p in [anchor] + positions
            )
            alt = gene.genomic_base(anchor_index)
            indels.append((anchor, ref, alt))
            truth.append(dict(variant_id=f"{CHROM}:{anchor}:{ref}>{alt}",
                              gene=name, consequence="inframe_deletion"))

    # Input that has to be rejected rather than annotated.
    gene1 = genes["GENE1"]
    bad_pos = gene1.positions[90]
    real_base = gene1.genomic_base(90)
    wrong_ref = next(b for b in "ACGT" if b != real_base)
    wrong_alt = next(b for b in "ACGT" if b not in (real_base, wrong_ref))
    snvs.append((bad_pos, wrong_ref, wrong_alt))
    rejected = [dict(variant_id=f"{CHROM}:{bad_pos}:{wrong_ref}>{wrong_alt}",
                     reason="reference allele does not match the reference sequence")]

    intron_pos = gene1.blocks[0][1] + 100
    snvs.append((intron_pos, "A", "G"))
    rejected.append(dict(variant_id=f"{CHROM}:{intron_pos}:A>G",
                         reason="not annotated on any transcript"))

    edge = gene1.blocks[0][1]  # last base of the first exon
    edge_ref = gene1.genomic_base(gene1.positions.index(edge))
    indels.append((edge, edge_ref + "AAA", edge_ref))
    # The pipeline's own wording for this, so the checker can look it up.
    rejected.append(dict(variant_id=f"{CHROM}:{edge}:{edge_ref}AAA>{edge_ref}",
                         reason="outside the coding blocks"))

    return snvs, indels, truth, rejected


def write_vcf(path, records):
    with gzip.open(path, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("##reference=file:/toy_GRCh38.fa\n")
        handle.write('##INFO=<ID=TVAF,Number=1,Type=Float,Description="tumour VAF">\n')
        handle.write(f"##contig=<ID={CHROM}>\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for pos, ref, alt in sorted(records):
            handle.write(
                f"{CHROM}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\tSOMATIC;TVAF=0.500\n")


def main(out_dir):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(SEED)

    genes, planted = build_genes(rng)
    snvs, indels, truth, rejected = make_variants(genes, planted)

    with gzip.open(out / "toy.gtf.gz", "wt") as handle:
        handle.write("#!toy annotation for testing\n")
        for gene in genes.values():
            handle.writelines(gene.gtf_lines())

    with gzip.open(out / "toy.cds.fa.gz", "wt") as handle:
        for gene in genes.values():
            handle.write(f">{gene.transcript} cds chromosome:TEST:{CHROM} "
                         f"gene:{gene.name}.1 gene_symbol:{gene.name}\n{gene.cds}\n")

    with gzip.open(out / "toy.pep.fa.gz", "wt") as handle:
        for gene in genes.values():
            handle.write(f">{gene.protein_id} pep chromosome:TEST:{CHROM} "
                         f"gene:{gene.name}.1 transcript:{gene.transcript} "
                         f"gene_symbol:{gene.name}\n{gene.protein}\n")

    write_vcf(out / "toy_snv.vcf.gz", snvs)
    write_vcf(out / "toy_indel.vcf.gz", indels)

    with open(out / "expression.tsv", "w") as handle:
        handle.write("gene_id\tgene_name\ttpm\n")
        for gene in genes.values():
            handle.write(f"{gene.name}\t{gene.name}\t{gene.tpm:g}\n")

    # What the pipeline is expected to produce.
    gene1 = genes["GENE1"]
    truth_file = dict(
        chromosome=CHROM,
        variants=truth,
        rejected=rejected,
        no_peptides_from=[
            entry["variant_id"] for entry in truth
            if entry["consequence"] in ("synonymous", "stop_gained")
        ],
        self_filtered_peptide=gene1.protein[60:69],
        gene_not_expressed="GENE4",
        positive_controls=["NLVPMVATV", "GILGFVFTL"],
    )
    (out / "truth.json").write_text(json.dumps(truth_file, indent=2) + "\n")

    print(f"Wrote {len(genes)} toy genes, {len(snvs)} SNVs and {len(indels)} indels to {out}/")
    print(f"  {len(truth)} variants with an expected consequence")
    print(f"  {len(rejected)} variants that have to be rejected")
    print(f"  self-filter decoy peptide: {truth_file['self_filtered_peptide']}")


if __name__ == "__main__":
    main(sys.argv[1])
