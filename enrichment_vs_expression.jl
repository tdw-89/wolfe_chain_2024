include("prelude.jl")

using .EnrichmentUtils

using Serialization
using JSON
using GLM

# function for serializing enrichment vectors to JSON for use in R
function serialize_to_json(file_path, vecs)
    mark_names = ["K27ac", "K4me3", "K9me3", "ATAC"]
    
    open(file_path, "w") do file
        JSON.print(file, Dict([mark => vs for (mark, vs) in zip(mark_names, vecs)]), 4)

    end
end

universal_range = GeneRange(TSS(), TES(), -500, 500)

# Peak files
chip_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_1_ensembl52/"
atac_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_2_ensembl52/"

# dicty ensembl 52 genome data
gff_data = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"

# Expression data
expr_data_file = "../../dicty_data/filtered/expr_data_filt_kallisto_ensembl52_single.tsv"

# Save directory
ser_data_dir = "../../dicty_data/julia_serialized/"

# MAIN SCRIPT: #

# Load genome data
ref_genome = loadgenome(gff_data, chrom_lengths_file)

# Load peak data
peak_files = [filter(fn -> !contains(fn, r"_S[AB]_"), readdir(chip_peak_file_dir, join=true)); readdir(atac_peak_file_dir, join=true)]
filter!(fn -> endswith(fn, ".narrowPeak"), peak_files)
peak_data = binpeaks(peak_files, chrom_lengths_file)
addtogenes!(ref_genome, peak_data)

# Load expression data
expr_data = CSV.read(expr_data_file, DataFrame)

# Rename the expression columns to a letter indicating with life-cycle stage the sample
# was taken from.
rename!(expr_data, ["GeneID", "V", "S", "M", "F"])

# Remove 'S' stage expression data, and rearrange so that the samples are in the order they
# will be for the peak data
# select!(expr_data, ["GeneID", "F", "M", "V"])
# Average gene expression counts across all life-cycle stages
insertcols!(expr_data, :Avg => mean.(eachrow(expr_data[:, 2:end])))
select!(expr_data, ["GeneID", "Avg"])

# log-transform data with an added pseudocount of 0.5
expr_data.Avg = log.(expr_data.Avg .+ 0.5)

# Add expression data to genome object
addexpression!(ref_genome, expr_data)

# Filter to genes that have expression data and have 500 bp upstream and 500 downstream before a chromosome end
filtered_gene_list = [gene for gene in ref_genome.genes[2] if has_expr(gene) && siginrange(gene, GeneRange(TSS(), TES(), -500, 500), peak_data=true)]
open("../../dicty_data/filtered/final_gene_list.txt", "w") do file
    for gene in filtered_gene_list
        println(file, gene.id)
    end
end

# Plot enrichment vs expression
k27ac_inds = [1, 4, 7]
k4me3_inds = [2, 5, 8]
k9me3_inds = [3, 6, 9]
atac_inds = [10, 11, 12]
sample_inds_vec = [k27ac_inds, k4me3_inds, k9me3_inds, atac_inds]

# Global signal distributions (mean/std-dev) for each sample group across all filtered genes
get_mean_sig(gene, gene_range, sample_inds) = mean([mean(getsiginrange(gene, gene_range, sample_ind)) for sample_ind in sample_inds])
get_global_mean(gene_list, gene_range, sample_inds) = mean([get_mean_sig(gene, gene_range, sample_inds) for gene in gene_list])
get_global_std_dev(gene_list, gene_range, sample_inds) = std([get_mean_sig(gene, gene_range, sample_inds) for gene in gene_list])
get_sig_dist(gene_list, gene_range, sample_inds) =
    EnrichmentUtils.SignalDistribution(
        get_global_mean(gene_list, gene_range, sample_inds),
        get_global_std_dev(gene_list, gene_range, sample_inds)
    )

global_means_vec = [get_sig_dist(filtered_gene_list, universal_range, sample_inds) for sample_inds in sample_inds_vec]

tss_enrich = plot_enrich_expr_region(expr_data, 
    filtered_gene_list, 
    sample_inds_vec, 
    [GeneRange(TSS(), TSS(), -500, 0) for i in 1:5], 
    fold_change_over_mean=false, 
    global_means=global_means_vec,
    z_min=z_min,
    z_max=z_max,
    return_figs=true)

serialize(joinpath(ser_data_dir, "tss_enrich_plots_expr.jls"), tss_enrich)

body_enrich = plot_enrich_expr_percent(expr_data, 
    filtered_gene_list, 
    sample_inds_vec, 
    fold_change_over_mean=false, 
    global_means=global_means_vec,
    z_min=z_min,
    z_max=z_max,
    return_figs=true)

serialize(joinpath(ser_data_dir, "body_enrich_plots_expr.jls"), body_enrich)

tes_enrich = plot_enrich_expr_region(expr_data,
    filtered_gene_list,
    sample_inds_vec,
    [GeneRange(TES(), TES(), 0, 500) for i in 1:5],
    fold_change_over_mean=false,
    global_means=global_means_vec,
    z_min=z_min,
    z_max=z_max,
    return_figs=true)

serialize(joinpath(ser_data_dir, "tes_enrich_plots_expr.jls"), tes_enrich)

bar_plots, kw_tests, means_vecs = plot_bar_expr(
    expr_data, 
    filtered_gene_list, 
    sample_inds_vec, 
    [universal_range, universal_range, universal_range, universal_range],
    global_means_vec, 
    [z_min, z_max];
    fold_change_over_mean=false,
    return_figs=true,
    horizontal=false
    )
eta_sq_vals = η².(kw_tests)
p_vals_perm_cor = [
    get_cor_expr(
        expr_data,
        universal_range,
        sample_ind,
        ref_genome
    ) for sample_ind in sample_inds_vec
]
adj_pvals_kw = adjust(pvalue.(kw_tests), BenjaminiHochberg())
adj_pvals_cor = adjust([pair[1] for pair in p_vals_perm_cor], BenjaminiHochberg())
serialize(joinpath(ser_data_dir, "bar_plots_expr.jls"), bar_plots)
serialize_to_json(joinpath(ser_data_dir, "means_vecs_expr.json"), means_vecs)