#!/usr/bin/env bash
# ============================================================================
# extract_gene2go.sh
# Extract organism-specific GO annotations from NCBI's gene2go file
#
# Usage:
#   bash extract_gene2go.sh <TAXID> [input_file] [output_file]
#
# Examples:
#   bash extract_gene2go.sh 310915                          # P. hypophthalmus
#   bash extract_gene2go.sh 8128                            # Nile tilapia
#   bash extract_gene2go.sh 310915 gene2go my_output.tsv   # custom paths
#
# Common aquaculture TaxIDs:
#   310915  Pangasianodon hypophthalmus (striped catfish)
#   7998    Ictalurus punctatus (channel catfish)
#   8128    Oreochromis niloticus (Nile tilapia)
#   29159   Crassostrea gigas (Pacific oyster)
#   8030    Salmo salar (Atlantic salmon)
#   8049    Oncorhynchus mykiss (rainbow trout)
# ============================================================================

set -euo pipefail

TAXID="${1:?Usage: bash extract_gene2go.sh <TAXID> [input_file] [output_file]}"
INPUT="${2:-gene2go}"
OUTPUT="${3:-gene2go_${TAXID}.tsv}"

# ── Check input exists ──
if [ ! -f "$INPUT" ] && [ ! -f "${INPUT}.gz" ]; then
    echo "Input file not found: $INPUT"
    echo "Download it first:"
    echo "  wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz"
    echo "  gunzip gene2go.gz"
    exit 1
fi

# ── Decompress on-the-fly if only .gz exists ──
if [ ! -f "$INPUT" ] && [ -f "${INPUT}.gz" ]; then
    echo "Found ${INPUT}.gz — extracting on-the-fly (no disk decompression needed)"
    echo "Filtering TaxID ${TAXID} ..."
    # Print header + matching rows
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
TERMS=$(awk -F'\t' 'NR > 1 { print $3 }' "$OUTPUT" | sort -u | wc -l)

echo "Done: ${OUTPUT}"
echo "  ${LINES} lines (including header)"
echo "  ${GENES} unique genes"
echo "  ${TERMS} unique GO terms"
