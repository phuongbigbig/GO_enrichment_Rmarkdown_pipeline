#!/usr/bin/env bash
# ============================================================================
# extract_gene2go.sh
# Extract organism-specific GO annotations from NCBI's gene2go file
#
# Usage:
#   bash extract_gene2go.sh <TAXID> [input_file] [output_file]
#
# The input_file can be:
#   - gene2go        (uncompressed)
#   - gene2go.gz     (compressed — decompressed on-the-fly)
#   - omitted        (auto-detects gene2go or gene2go.gz in current dir)
# ============================================================================

set -euo pipefail

TAXID="${1:?Usage: bash extract_gene2go.sh <TAXID> [input_file] [output_file]}"
INPUT="${2:-gene2go}"
OUTPUT="${3:-gene2go_${TAXID}.tsv}"

# ── Resolve input file ──
if [ -f "$INPUT" ]; then
    RESOLVED="$INPUT"
elif [ -f "${INPUT}.gz" ]; then
    RESOLVED="${INPUT}.gz"
else
    echo "Error: input file not found: $INPUT (also tried ${INPUT}.gz)"
    echo "Download it first:"
    echo "  wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
    exit 1
fi

echo "Filtering ${RESOLVED} for TaxID ${TAXID} ..."

# ── Use gunzip -c for .gz (works on both macOS and Linux), plain cat otherwise ──
if [[ "$RESOLVED" == *.gz ]]; then
    gunzip -c "$RESOLVED"
else
    cat "$RESOLVED"
fi | LC_ALL=C awk -F'\t' -v tid="$TAXID" '
    NR == 1 { gsub(/^#/, "", $0); print; next }
    /^#/    { next }
    $1 == tid
' > "$OUTPUT"

# ── Report ──
LINES=$(wc -l < "$OUTPUT")
GENES=$(LC_ALL=C awk -F'\t' 'NR > 1 { print $2 }' "$OUTPUT" | sort -u | wc -l)
TERMS=$(LC_ALL=C awk -F'\t' 'NR > 1 { print $3 }' "$OUTPUT" | sort -u | wc -l)

echo "Done: ${OUTPUT}"
echo "  ${LINES} lines (including header)"
echo "  ${GENES} unique genes"
echo "  ${TERMS} unique GO terms"
