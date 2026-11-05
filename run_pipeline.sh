#!/usr/bin/env bash
# Run the neoantigen pipeline.
#
#   ./run_pipeline.sh --test   simulated data, then check the results -> results/test/
#   ./run_pipeline.sh          HCC1395 (SEQC2 truth set)              -> results/
#
set -euo pipefail

SETTINGS=config/settings.yaml
LENGTHS=$(python - "$SETTINGS" <<'PY'
import sys, yaml
print(",".join(str(x) for x in yaml.safe_load(open(sys.argv[1]))["peptide_lengths"]))
PY
)

if [[ "${1:-}" == "--test" ]]; then
    MODE=test
    DATA=data/test
    REF=data/test/ref
    OUT=results/test
    HLA=config/hla_test.tsv
else
    MODE=real
    DATA=data/raw
    REF=data/ref
    OUT=results
    HLA=config/hla.tsv
fi
mkdir -p "$DATA" "$REF" "$OUT"

# The MHCflurry models are a separate download of about 100 MB.
if ! mhcflurry-downloads path models_class1_presentation > /dev/null 2>&1; then
    echo "== Fetching the MHCflurry models (once)"
    mhcflurry-downloads fetch models_class1_presentation
fi

if [[ "$MODE" == test ]]; then
    echo "== 1. Simulate data with a known answer"
    python scripts/simulate_test_data.py "$DATA"
    GTF="$DATA/toy.gtf.gz"; CDS="$DATA/toy.cds.fa.gz"; PEP="$DATA/toy.pep.fa.gz"
    VCFS=("$DATA/toy_snv.vcf.gz" "$DATA/toy_indel.vcf.gz")
else
    echo "== 1. Download the reference, the variants and the expression table"
    python scripts/01_fetch_data.py "$SETTINGS" "$DATA"
    GTF="$DATA/Homo_sapiens.GRCh38.110.chr.gtf.gz"
    CDS="$DATA/Homo_sapiens.GRCh38.cds.all.fa.gz"
    PEP="$DATA/Homo_sapiens.GRCh38.pep.all.fa.gz"
    VCFS=("$DATA"/high-confidence_sSNV_in_HC_regions_*.vcf.gz
          "$DATA"/high-confidence_sINDEL_in_HC_regions_*.vcf.gz)
fi

echo "== 2. Build the transcript database"
python scripts/02_build_reference.py "$GTF" "$CDS" "$PEP" "$REF"

echo "== 3. Annotate the variants and translate the mutant proteins"
python scripts/03_annotate_variants.py "$REF" "$OUT" "${VCFS[@]}"

echo "== 4. Cut out the mutant peptides and remove the ones that are not tumour-specific"
python scripts/04_peptides.py "$REF" "$OUT" "$LENGTHS"

echo "== 5. Predict presentation by the sample's HLA alleles"
python scripts/05_predict_binding.py "$HLA" "$OUT"

echo "== 6. Filter, rank and draw the figures"
python scripts/06_rank.py "$SETTINGS" "$DATA/expression.tsv" "$OUT"

if [[ "$MODE" == test ]]; then
    echo "== 7. Check the results against the planted truth"
    python scripts/check_test_results.py "$DATA" "$OUT"
fi
echo "Done. Results are in $OUT/"
