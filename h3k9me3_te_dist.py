import os
import subprocess
import sys

BEDTOOLS_DIR = "/Users/thomaswolfe/Documents/bedtools2/bin"
os.environ["PATH"] = f"{BEDTOOLS_DIR}:{os.environ.get('PATH', '')}"
genes_bed = "../../dicty_data/te_density/cds_coordinates_sorted_3_column.bed"
tes_bed = "../../dicty_data/te_density/te_coordinates_sorted_merged.bed"
peaks_bed = "../../dicty_data/te_density/h3k9me3_sorted_merged.bed"
gene_te_coords_only = "../../gene_te_coords_only.bed"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
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

    print(f"Total TE bases: {total_te_bases:,}")
    print(f"Overlapping bases: {overlap_te_bases:,}")
    print(f"Percentage of TE bases with H3K9me3: {pct_overlap_te_by_k9:.2f}%\n")
    print(f"Percentage of H3K9me3 bases overlapping a TE: {pct_overlap_k9_by_te:.2f}%\n")
    print(f"Total non-TE gene bases: {total_gene_bases:,}")
    print(f"Overlapping bases: {overlap_gene_bases:,}")
    print(f"Percentage of non-TE gene bases with H3K9me3: {pct_overlap_gene_by_k9:.2f}%\n")
    print(f"Percentage of H3K9me3 bases overlapping a non-TE gene: {pct_overlap_k9_by_gene:.2f}%\n")

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
            f"bedtools shuffle -i {tes_bed} -g {chrom_lengths_file} -incl {gene_te_coords_only} | "
            f"bedtools intersect -a stdin -b {peaks_bed} | "
            f"bedtools sort -i - | bedtools merge -i - | "
            f"awk '{{sum += $3 - $2}} END {{print sum+0}}'"
        )
        
        rand_overlap = int(run_command(shuffle_cmd))
        
        if rand_overlap >= obs_overlap:
            count_greater_or_equal += 1
            
        if i % 100 == 0:
            print(f"   ...completed {i}/{iterations} permutations")

    print("\n4. Calculating empirical p-value...")
    p_value = count_greater_or_equal / iterations

    print("-" * 40)
    print("RESULTS:")
    print(f"Total iterations: {iterations}")
    print(f"Random iterations >= Observed: {count_greater_or_equal}")
    print(f"Empirical P-value: {p_value}")
    
    if p_value == 0.0:
        print("Note: Report as p < 0.001 (the random overlap never reached the observed overlap).")
    print("-" * 40)

if __name__ == "__main__":
    main()