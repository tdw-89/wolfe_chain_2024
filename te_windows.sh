#!/usr/bin/env bash
# Build the BED inputs for te_density.jl and h3k9me3_te_dist.py from the files written by te_distance.jl,
# te_distance_human.jl and h3k9me3_dist_dicty.jl: sort the genes, sort and merge the TEs, then list every
# TE within 100 kb of each gene.
#
# Run from anywhere after those scripts; run_analysis.sh calls it in order.
set -euo pipefail

BEDTOOLS="${BEDTOOLS:-$(command -v bedtools || echo "$HOME/Documents/bedtools2/bin/bedtools")}"
WINDOW=100000   # must match WINDOW in te_density.jl

cd "$(dirname "$0")/../../dicty_data/te_density"
for tag in "" "_human"; do
    "$BEDTOOLS" sort -i "cds_coordinates${tag}.bed" > "cds_coordinates${tag}_sorted.bed"
    "$BEDTOOLS" sort -i "te_coordinates${tag}.bed" > "te_coordinates${tag}_sorted.bed"
    "$BEDTOOLS" merge -i "te_coordinates${tag}_sorted.bed" > "te_coordinates${tag}_sorted_merged.bed"
    "$BEDTOOLS" window -w "$WINDOW" -a "cds_coordinates${tag}_sorted.bed" -b "te_coordinates${tag}_sorted_merged.bed" \
        > "tes_in_window${tag}.bed"
done

# D. discoideum genes and H3K9me3 peaks for h3k9me3_te_dist.py
cut -f1-3 cds_coordinates_sorted.bed > cds_coordinates_sorted_3_column.bed
"$BEDTOOLS" sort -i h3k9me3_union.narrowPeak | "$BEDTOOLS" merge -i - > h3k9me3_sorted_merged.bed
