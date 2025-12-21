include("prelude.jl")

using Serialization
using JSON
using Random

# function for serializing enrichment vectors to JSON for use in R
function serialize_to_json(file_path, vecs)
    mark_names = ["K27ac", "K4me3", "K9me3", "ATAC"]
    
    open(file_path, "w") do file
        JSON.print(file, Dict([mark => vs for (mark, vs) in zip(mark_names, vecs)]), 4)

    end
end

# Peak files
chip_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_1_ensembl52/"
atac_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_2_ensembl52/"
sig_region_file = "../../dicty_data/sig_regions.csv"

# Genome data
gff_source = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
paralog_file = "../../dicty_data/filtered/paralog_filt.tsv"
ensembl_cds_id_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/cds_ids.txt"
blacklist_file = "../../dicty_data/blacklist_with_tes.csv"
final_gene_list = "../../dicty_data/filtered/final_gene_list.txt"
singleton_file = "../../dicty_data/filtered/singleton_filt.tsv"

# Save directory
ser_data_dir = "../../dicty_data/julia_serialized/"

# MAIN SCRIPT: #

# Load genome data
ref_genome = loadgenome(gff_source, chrom_lengths_file)

# Load final filtered gene list
final_gene_list = open(final_gene_list) do file
    readlines(file)
end

# Load peak data
peak_files = [readdir(chip_peak_file_dir, join=true); readdir(atac_peak_file_dir, join=true)]
peak_files = filter(fn -> endswith(fn, ".narrowPeak"), peak_files)
peak_files = filter(fn -> !contains(fn, r"_S[AB]+_"), peak_files)
peak_data = binpeaks(peak_files, chrom_lengths_file)
addtogenes!(ref_genome, peak_data)

# Load paralog data
paralog_data = CSV.read(paralog_file, DataFrame)
filter!(row -> row["dS"] <= 3, paralog_data)
CSV.write("../../dicty_data/filtered/paralog_ds_filt.csv", paralog_data)
select!(paralog_data, ["GeneID", "ParalogID", "dS"])

# paralog_data[findall(map(id -> !siginrange(get(ref_genome, id),  GeneRange(TSS(), TES(), -500, 500)), paralog_data.ParalogID)),:]

# select!(paralog_data, ["GeneID", "ParalogID", "MaxPerc"])
paralog_ids = vcat(paralog_data[:, 1], paralog_data[:, 2])

# Filter to genes that have 500 bp upstream and 250 downstream before a chromosome end
filtered_gene_list = [gene for gene in ref_genome.genes[2] if gene.id in final_gene_list]
filtered_paralog_list = [gene for gene in filtered_gene_list if gene.id in paralog_ids]

# Plot enrichment vs expression
k27ac_inds = [1, 4, 7]
k4me3_inds = [2, 5, 8]
k9me3_inds = [3, 6, 9]
atac_inds = [10, 11, 12]
sample_inds_vec = [k27ac_inds, 
                   k4me3_inds, 
                   k9me3_inds, 
                   atac_inds]
# get the global mean enrichment for each sample across all filtered genes
global_means_vec = 
[mean([mean([mean(getsiginrange(gene, GeneRange(TSS(), TES(), -500, 500), sample_ind)) for sample_ind in sample_inds]) for gene in filtered_gene_list]) for sample_inds in sample_inds_vec]

tss_enrich = plot_enrich_region(
    paralog_data, 
    filtered_paralog_list, 
    sample_inds_vec, 
    [GeneRange(TSS(), TSS(), -500, 0) for i in 1:4], 
    fold_change_over_mean=true, 
    global_means=global_means_vec,
    z_min=0,
    z_max=4,
    return_figs=true
    )

serialize(joinpath(ser_data_dir, "tss_enrich_plots_dS.jls"), tss_enrich)

body_enrich = plot_enrich_percent(paralog_data, 
filtered_paralog_list, 
    sample_inds_vec, 
    fold_change_over_mean=true, 
    global_means=global_means_vec,
    z_min=0,
    z_max=4,
    return_figs=true)

serialize(joinpath(ser_data_dir, "body_enrich_plots_dS.jls"), body_enrich)

tes_enrich = plot_enrich_region(paralog_data,
filtered_paralog_list,
    sample_inds_vec,
    [GeneRange(TES(), TES(), 0, 500) for i in 1:4],
    fold_change_over_mean=true,
    global_means=global_means_vec,
    z_min=0,
    z_max=4,
    return_figs=true)

serialize(joinpath(ser_data_dir, "tes_enrich_plots_dS.jls"), tes_enrich)

# Plot bar plots of coverage in significant regions
sig_region_df = CSV.read(sig_region_file, DataFrame)
bar_plots, kw_tests, means_vecs = plot_bar(paralog_data, filtered_paralog_list, sample_inds_vec, [GeneRange(TSS(), TES(), sig_region_df.Start[1], sig_region_df.End[1]), # K27ac
                                                                            GeneRange(TSS(), TES(), sig_region_df.Start[2], sig_region_df.End[2]), # K4me3
                                                                            GeneRange(TSS(), TES(), sig_region_df.Start[3], sig_region_df.End[3]), # K9me3
                                                                            GeneRange(TSS(), TES(), sig_region_df.Start[4], sig_region_df.End[4])], # ATAC
                                                            global_means_vec, [0,4], true, true);
p_vals_perm_cor = [get_cor(paralog_data, gene_range, sample_ind, ref_genome, global_mean) for (sample_ind, gene_range, global_mean) in zip(sample_inds_vec,
                            [GeneRange(TSS(), TES(), sig_region_df.Start[1], sig_region_df.End[1]),
                             GeneRange(TSS(), TES(), sig_region_df.Start[2], sig_region_df.End[2]),
                             GeneRange(TSS(), TES(), sig_region_df.Start[3], sig_region_df.End[3]),
                             GeneRange(TSS(), TES(), sig_region_df.Start[4], sig_region_df.End[4])],
                             global_means_vec)]
adj_p_vals_perm_cor = adjust([pair[1] for pair in p_vals_perm_cor], Bonferroni())
adj_p_vals_kw = adjust(pvalue.(kw_tests), Bonferroni())
serialize(joinpath(ser_data_dir, "bar_plots_dS.jls"), bar_plots)
serialize_to_json(joinpath(ser_data_dir, "means_vecs_ds.json"), means_vecs)