#!/usr/bin/env bash
# ============================================================================
# extract_gene_info.sh
# Extract organism-specific gene info from NCBI's gene_info file
# (used for Symbol ↔ GeneID mapping in the GO enrichment pipeline)
#
# Usage:
#   bash extract_gene_info.sh <TAXID> [input_file] [output_file]
#
# Examples:
#   bash extract_gene_info.sh 310915
#   bash extract_gene_info.sh 310915 gene_info my_gene_info.tsv
# ============================================================================

set -euo pipefail

TAXID="${1:?Usage: bash extract_gene_info.sh <TAXID> [input_file] [output_file]}"
INPUT="${2:-gene_info}"
OUTPUT="${3:-gene_info_${TAXID}.tsv}"

# ── Check input exists ──
if [ ! -f "$INPUT" ] && [ ! -f "${INPUT}.gz" ]; then
    echo "Input file not found: $INPUT"
    echo "Download it first:"
    echo "  wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz"
    echo "  gunzip gene_info.gz"
    exit 1
fi

# ── Decompress on-the-fly if only .gz exists ──
if [ ! -f "$INPUT" ] && [ -f "${INPUT}.gz" ]; then
    echo "Found ${INPUT}.gz — extracting on-the-fly"
    echo "Filtering TaxID ${TAXID} ..."
    zcat "${INPUT}.gz" | awk -F'\t' -v tid="$TAXID" '
        NR == 1 { gsub(/^#/, "", $0); print; next }
        /^#/    { next }
        $1 == tid
    ' > "$OUTPUT"
else
    echo "Filtering ${INPUT} for TaxID ${TAXID} ..."
    awk -F'\t' -v tid="$TAXID" '
        NR == 1 { gsub(/^#/, "", $0); print; next }
        /^#/    { next }
        $1 == tid
    ' "$INPUT" > "$OUTPUT"
fi

# ── Report ──
LINES=$(wc -l < "$OUTPUT")
GENES=$(awk -F'\t' 'NR > 1 { print $2 }' "$OUTPUT" | sort -u | wc -l)

echo "Done: ${OUTPUT}"
echo "  ${LINES} lines (including header)"
echo "  ${GENES} unique genes"
