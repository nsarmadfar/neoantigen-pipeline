"""Small helpers shared by the pipeline steps: the genetic code and FASTA."""
import gzip

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# The 20 amino acids. Anything else (X from an N in the sequence, or U from
# selenocysteine) means the peptide cannot be scored, so those are dropped.
AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")


def translate(seq, to_stop=True):
    """Translate a coding sequence.

    Returns the protein and whether a stop codon was reached. Codons that
    contain an N become X, so they can be filtered out later instead of
    quietly becoming a wrong residue.
    """
    protein = []
    for i in range(0, len(seq) - len(seq) % 3, 3):
        residue = CODON_TABLE.get(seq[i:i + 3].upper(), "X")
        if residue == "*" and to_stop:
            return "".join(protein), True
        protein.append(residue)
    return "".join(protein), False


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def strip_version(identifier):
    """ENST00000361390.2 -> ENST00000361390"""
    return identifier.split(".")[0]


def read_fasta(path):
    """Yield (header, sequence) from a plain or gzipped FASTA file."""
    opener = gzip.open if str(path).endswith(".gz") else open
    header, chunks = None, []
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].rstrip("\n"), []
            else:
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def read_table(path):
    """Read a tab-separated file with a header row into a list of dicts."""
    rows = []
    with open(path) as handle:
        columns = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            if not line.strip():
                continue
            rows.append(dict(zip(columns, line.rstrip("\n").split("\t"))))
    return rows


def write_table(path, rows, columns=None):
    """Write a list of dicts as a tab-separated file."""
    columns = columns or (list(rows[0].keys()) if rows else [])
    with open(path, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(c, "")) for c in columns) + "\n")
