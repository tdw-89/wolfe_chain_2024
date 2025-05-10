using Serialization
using MultipleTesting
using JSON

# Custom lib src:
include("./custom_lib/load_gff.jl")
include("./custom_lib/genomic_data.jl")
include("./custom_lib/enrichment_utils.jl")

# function for serializing enrichment vectors to JSON for use in R
function serialize_to_json(file_path, vecs)
    mark_names = ["K27ac", "K4me3", "K9me3", "ATAC"]
    
    open(file_path, "w") do file
        JSON.print(file, Dict([mark => vs for (mark, vs) in zip(mark_names, vecs)]), 4)

    end
end

# Peak files
chip_peak_file_dir = "../../../../data/wang_et_al/processed/run_1_ensembl52/"
atac_peak_file_dir = "../../../../data/wang_et_al/processed/run_2_ensembl52/"
sig_region_file = "./data/sig_regions.csv"

# dicty ensembl 52 genome data
gff_data = "../../../../data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../../../data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"

# Expression data
expr_data_file = "./data/filtered/expr_data_filt_kallisto_ensembl52_single.tsv"

# Save directory
ser_data_dir = "./data/julia_serialized/"

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
open("./data/filtered/final_gene_list.txt", "w") do file
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
global_means_vec = 
[mean([mean([mean(getsiginrange(gene, GeneRange(TSS(), TES(), -500, 500), sample_ind)) for sample_ind in sample_inds]) for gene in filtered_gene_list]) for sample_inds in sample_inds_vec]

tss_enrich = plot_enrich_expr_region(expr_data, 
    filtered_gene_list, 
    sample_inds_vec, 
    [GeneRange(TSS(), TSS(), -500, 0) for i in 1:5], 
    fold_change_over_mean=true, 
    global_means=global_means_vec,
    z_min=0,
    z_max=4,
    return_figs=true)

serialize(joinpath(ser_data_dir, "tss_enrich_plots_expr.jls"), tss_enrich)

body_enrich = plot_enrich_expr_percent(expr_data, 
    filtered_gene_list, 
    sample_inds_vec, 
    fold_change_over_mean=true, 
    global_means=global_means_vec,
    z_min=0,
    z_max=4,
    return_figs=true)

serialize(joinpath(ser_data_dir, "body_enrich_plots_expr.jls"), body_enrich)

tes_enrich = plot_enrich_expr_region(expr_data,
    filtered_gene_list,
    sample_inds_vec,
    [GeneRange(TES(), TES(), 0, 500) for i in 1:5],
    fold_change_over_mean=true,
    global_means=global_means_vec,
    z_min=0,
    z_max=4,
    return_figs=true)

serialize(joinpath(ser_data_dir, "tes_enrich_plots_expr.jls"), tes_enrich)

sig_region_df = CSV.read(sig_region_file, DataFrame)
bar_plots, kw_tests, means_vecs = plot_bar_expr(expr_data, filtered_gene_list, sample_inds_vec, [GeneRange(TSS(), TES(), sig_region_df.Start[1], parse(Int, sig_region_df.End[1])), # K27ac
                                                                            GeneRange(TSS(), TES(), sig_region_df.Start[2], parse(Int, sig_region_df.End[2])), # K4me3
                                                                            GeneRange(TSS(), TSS(), sig_region_df.Start[3], 0), # K9me3
                                                                            GeneRange(TSS(), TES(), sig_region_df.Start[4], parse(Int, sig_region_df.End[4]))], # ATAC
                                                        global_means_vec, [0,4], true, true);
p_vals_perm_cor = [get_cor_expr(expr_data, 
                                         gene_range, 
                                         sample_ind, 
                                         ref_genome) for (sample_ind, gene_range) in zip(sample_inds_vec,
                            [GeneRange(TSS(), TES(), sig_region_df.Start[1], parse(Int, sig_region_df.End[1])),
                             GeneRange(TSS(), TES(), sig_region_df.Start[2], parse(Int, sig_region_df.End[2])),
                             GeneRange(TSS(), TES(), 0, 0),
                             GeneRange(TSS(), TES(), sig_region_df.Start[4], parse(Int, sig_region_df.End[4]))])]
adj_pvals_kw = adjust(pvalue.(kw_tests), Bonferroni())
adj_pvals_cor = adjust([pair[1] for pair in p_vals_perm_cor], Bonferroni())
serialize(joinpath(ser_data_dir, "bar_plots_expr.jls"), bar_plots)
serialize_to_json(joinpath(ser_data_dir, "means_vecs_expr.json"), means_vecs)