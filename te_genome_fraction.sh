#!/usr/bin/env bash
# TE content of each genome (Table 1): non-redundant TE bases by class as a percentage of the assembly
# length. Class II = DNA transposons, Class I = all other classified TEs. D. discoideum TE genes from the
# annotation (type "_converted_") have no class, so they count only towards "All (with TE genes)".
#
# Run after te_distance.jl and te_distance_human.jl. Writes ../../dicty_data/stats/te_genome_fraction.csv.
set -euo pipefail

BEDTOOLS="${BEDTOOLS:-$(command -v bedtools || echo "$HOME/Documents/bedtools2/bin/bedtools")}"

cd "$(dirname "$0")/../../dicty_data"
mkdir -p stats
out=stats/te_genome_fraction.csv
echo "Species,Set,Bases,GenomeBases,Percent" > "$out"

# Write one row: species, set name, TE BED file, genome lengths file, awk filter on the TE type column
row() {
    local bases genome
    bases=$(awk -F'\t' "$5" "$3" | "$BEDTOOLS" sort -i - | "$BEDTOOLS" merge -i - | awk '{s += $3 - $2} END {print s + 0}')
    genome=$(awk '{s += $2} END {print s}' "$4")
    awk -v sp="$1" -v set="$2" -v b="$bases" -v g="$genome" \
        'BEGIN {printf "%s,%s,%d,%d,%.2f\n", sp, set, b, g, 100 * b / g}' >> "$out"
}

dicty_tes=te_density/te_coordinates.bed
dicty_genome=AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt
row "D. discoideum" "Class I" "$dicty_tes" "$dicty_genome" '$4 !~ /^DNA/ && $4 != "_converted_"'
row "D. discoideum" "Class II" "$dicty_tes" "$dicty_genome" '$4 ~ /^DNA/'
row "D. discoideum" "All (predicted TEs)" "$dicty_tes" "$dicty_genome" '$4 != "_converted_"'
row "D. discoideum" "All (with TE genes)" "$dicty_tes" "$dicty_genome" '1'

human_tes=te_density/te_coordinates_human.bed
human_genome=mammals/primates/h_sapiens/Ensembl_99/chromosome_lengths.txt
row "H. sapiens" "Class I" "$human_tes" "$human_genome" '$4 != "DNA"'
row "H. sapiens" "Class II" "$human_tes" "$human_genome" '$4 == "DNA"'
row "H. sapiens" "All" "$human_tes" "$human_genome" '1'
