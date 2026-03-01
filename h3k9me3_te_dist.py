import os
import subprocess
import sys


BEDTOOLS_DIR = "/Users/thomaswolfe/Documents/bedtools2/bin"
os.environ["PATH"] = f"{BEDTOOLS_DIR}:{os.environ.get('PATH', '')}"
genes_bed = "../../dicty_data/non_te_coding_genes.bed"
tes_bed = "../../dicty_data/predicted_tes.bed"
peaks_bed = "../../dicty_data/h3k9me3_peaks.bed"
gene_te_coords_only = "../../gene_te_coords_only.bed"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
iterations = 10000

def run_command(cmd):
    """
    Executes a shell command using subprocess, captures the output, 
    and handles errors gracefully.
    """
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
    
    print("Calculating bases overlapping H3K9me3 peaks...")
    overlap_bases_cmd = f"bedtools intersect -a {tes_bed} -b {peaks_bed} | bedtools sort -i - | bedtools merge -i - | awk '{{sum += $3 - $2}} END {{print sum+0}}'"
    overlap_bases = int(run_command(overlap_bases_cmd))
    
    if total_te_bases > 0:
        pct_overlap = (overlap_bases / total_te_bases) * 100
    else:
        pct_overlap = 0.0

    print(f"Total TE bases: {total_te_bases:,}")
    print(f"Overlapping bases: {overlap_bases:,}")
    print(f"Percentage of TE bases with H3K9me3: {pct_overlap:.2f}%\n")

    print("1. Creating the trimmed coordinates file for shuffling...")
    # Concatenate and sort to create the inclusion file for shuffling
    run_command(f"cat {genes_bed} {tes_bed} | bedtools sort > {gene_te_coords_only}")

    print("2. Calculating OBSERVED overlaps...")
    obs_cmd = f"bedtools intersect -a {tes_bed} -b {peaks_bed} -u | wc -l"
    obs_overlap = int(run_command(obs_cmd))
    print(f"   Observed TE overlaps: {obs_overlap}\n")

    print(f"3. Running {iterations} permutations")
    count_greater_or_equal = 0

    for i in range(1, iterations + 1):
        shuffle_cmd = (
            f"bedtools shuffle -i {tes_bed} -g {chrom_lengths_file} -incl {gene_te_coords_only} | "
            f"bedtools intersect -a stdin -b {peaks_bed} -u | wc -l"
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