"""Overlap of H3K9me3 peaks with TEs and non-TE coding genes in D. discoideum, and a permutation test of
whether TEs overlap H3K9me3 more than expected when shuffled within gene and TE regions.

Inputs are written by te_windows.sh. Results go to dicty_data/stats/h3k9me3_te_overlap.csv.
"""
import csv
import os
import subprocess
import sys

# bedtools location: $BEDTOOLS_DIR if set, otherwise ~/Documents/bedtools2/bin (also searched after PATH)
BEDTOOLS_DIR = os.environ.get("BEDTOOLS_DIR", os.path.expanduser("~/Documents/bedtools2/bin"))
os.environ["PATH"] = f"{BEDTOOLS_DIR}:{os.environ.get('PATH', '')}"
os.chdir(os.path.dirname(os.path.abspath(__file__)))
te_dir = "../../dicty_data/te_density"
genes_bed = f"{te_dir}/cds_coordinates_sorted_3_column.bed"
tes_bed = f"{te_dir}/te_coordinates_sorted_merged.bed"
peaks_bed = f"{te_dir}/h3k9me3_sorted_merged.bed"
gene_te_coords_only = f"{te_dir}/gene_te_coords_only.bed"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
stats_file = "../../dicty_data/stats/h3k9me3_te_overlap.csv"
iterations = 10000

def run_command(cmd):
    try:
        # shell=True allows us to use pipes (|) directly in the command string
        result = subprocess.run(cmd, shell=True, text=True, check=True, capture_output=True)
        return result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f"Error executing command:\n{cmd}")
        print(f"Error Message:\n{e.stderr}")
        sys.exit(1)

def main():
    print("Calculating total non-redundant TE bases...")
    # Note: Double curly braces {{ }} are used to escape the literal braces in Python f-strings
    total_te_cmd = f"bedtools sort -i {tes_bed} | bedtools merge -i - | awk '{{sum += $3 - $2}} END {{print sum+0}}'"
    total_te_bases = int(run_command(total_te_cmd))
    total_k9_cmd = f"bedtools sort -i {peaks_bed} | bedtools merge -i - | awk '{{sum += $3 - $2}} END {{print sum+0}}'"
    total_k9_bases = int(run_command(total_k9_cmd))
    total_genes_cmd = f"bedtools sort -i {genes_bed} | bedtools merge -i - | awk '{{sum += $3 - $2}} END {{print sum+0}}'"
    total_gene_bases = int(run_command(total_genes_cmd))

    print("Calculating TE bases overlapping H3K9me3 peaks...")
    overlap_te_bases_cmd = f"bedtools intersect -a {tes_bed} -b {peaks_bed} | bedtools sort -i - | bedtools merge -i - | awk '{{sum += $3 - $2}} END {{print sum+0}}'"
    overlap_te_bases = int(run_command(overlap_te_bases_cmd))
    print("Calculating non-TE gene bases overlapping H3K9me3 peaks...")
    overlap_gene_bases_cmd = f"bedtools intersect -a {genes_bed} -b {peaks_bed} | bedtools sort -i - | bedtools merge -i - | awk '{{sum += $3 - $2}} END {{print sum+0}}'"
    overlap_gene_bases = int(run_command(overlap_gene_bases_cmd))
    
    if total_te_bases > 0:
        pct_overlap_te_by_k9 = (overlap_te_bases / total_te_bases) * 100
        pct_overlap_k9_by_te = (overlap_te_bases / total_k9_bases) * 100
    else:
        pct_overlap_te_by_k9 = 0.0
        pct_overlap_k9_by_te = 0.0
    
    if total_gene_bases > 0:
        pct_overlap_gene_by_k9 = (overlap_gene_bases / total_gene_bases) * 100
        pct_overlap_k9_by_gene = (overlap_gene_bases / total_k9_bases) * 100
    else:
        pct_overlap_gene_by_k9 = 0.0
        pct_overlap_k9_by_gene = 0.0


    print("1. Creating the trimmed coordinates file for shuffling...")
    # Concatenate and sort to create the inclusion file for shuffling
    run_command(f"cat {genes_bed} {tes_bed} | bedtools sort > {gene_te_coords_only}")

    print("2. Calculating OBSERVED overlaps...")
    obs_overlap = overlap_te_bases
    print(f"   Observed TE overlap bases: {obs_overlap:,}\n")

    print(f"3. Running {iterations} permutations")
    count_greater_or_equal = 0

    for i in range(1, iterations + 1):
        shuffle_cmd = (
            f"bedtools shuffle -seed {i} -i {tes_bed} -g {chrom_lengths_file} -incl {gene_te_coords_only} | "
            f"bedtools intersect -a stdin -b {peaks_bed} | "
            f"bedtools sort -i - | bedtools merge -i - | "
            f"awk '{{sum += $3 - $2}} END {{print sum+0}}'"
        )
        
        rand_overlap = int(run_command(shuffle_cmd))
        
        if rand_overlap >= obs_overlap:
            count_greater_or_equal += 1
            
        if i % 100 == 0:
            print(f"   ...completed {i}/{iterations} permutations")

    # The observed arrangement counts as one permutation, so the smallest possible p is 1 / (iterations + 1)
    p_value = (count_greater_or_equal + 1) / (iterations + 1)

    os.makedirs(os.path.dirname(stats_file), exist_ok=True)
    with open(stats_file, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Statistic", "Value"])
        writer.writerows([
            ["TE bases", total_te_bases],
            ["H3K9me3 bases", total_k9_bases],
            ["Non-TE gene bases", total_gene_bases],
            ["TE bases overlapping H3K9me3", overlap_te_bases],
            ["Non-TE gene bases overlapping H3K9me3", overlap_gene_bases],
            ["Percent of TE bases with H3K9me3", pct_overlap_te_by_k9],
            ["Percent of H3K9me3 bases in TEs", pct_overlap_k9_by_te],
            ["Percent of non-TE gene bases with H3K9me3", pct_overlap_gene_by_k9],
            ["Percent of H3K9me3 bases in non-TE genes", pct_overlap_k9_by_gene],
            ["Permutations", iterations],
            ["Permutations with overlap >= observed", count_greater_or_equal],
            ["Permutation p", p_value],
        ])
    print(f"Results written to {stats_file}")

if __name__ == "__main__":
    main()