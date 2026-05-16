include("prelude.jl")

using .RepeatUtils
using .EnrichmentUtils

# Files:
gff_source = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
chip_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_1_ensembl52/"
te_id_file = "blacklists/ensembl_te_ids_with_predicted.tsv"
peak_files = filter(fn -> contains(fn, "K9me3") && !contains(fn, r"_S[AB]_"), readdir(chip_peak_file_dir, join=true))
peak_data = binpeaks(peak_files, chrom_lengths_file)

for i in 1:length(peak_data.samples[1].chroms)
    sig1 = peak_data.samples[1].chroms[i].signal
    sig2 = peak_data.samples[2].chroms[i].signal
    sig3 = peak_data.samples[3].chroms[i].signal
    sig_union = sig1 .| sig2 .| sig3
    peak_data.samples[1].chroms[i].signal = sig_union
end
sig_union = deepcopy(peak_data.samples[1])
peak_data = nothing

# Load genome data
ref_genome = loadgenome(gff_source, chrom_lengths_file)
te_ids = readlines(te_id_file)

peak_bed = DataFrame(:chrom=>String[], :start=>Int[], :end=>Int[])
for chrom in sig_union.chroms
    chrom_name = chrom.name
    signal = chrom.signal
    in_peak = false
    peak_start = 0
    for pos in 1:length(signal)
        if signal[pos] && !in_peak
            in_peak = true
            peak_start = pos
        elseif !signal[pos] && in_peak
            in_peak = false
            push!(peak_bed, (chrom_name, peak_start, pos - 1))
        end
    end
    # Handle case where we end in a peak
    if in_peak
        push!(peak_bed, (chrom_name, peak_start, length(signal)))
    end
end
CSV.write("../../dicty_data/", peak_bed, delim='\t', header=false)

# Create bed file for non-TE coding genes
gene_bed = DataFrame(:chrom=>String[], :start=>Int[], :end=>Int[])
for (chrom_name, chrom) in ref_genome.scaffolds
    for gene in chrom.genes
        if gene.id ∈ te_ids
            continue
        end
        gene_start = gene.gene_start
        gene_end = gene.gene_end
        push!(gene_bed, (chrom_name, gene_start, gene_end))
    end
end
CSV.write("../../dicty_data/non_te_coding_genes.bed", gene_bed, delim='\t', header=false)
