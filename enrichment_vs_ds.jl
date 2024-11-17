using Serialization
using JSON
using MultipleTesting

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
# chip_peak_file_dir = "../../../../data/wang_et_al/processed/run_1_ensembl52/"
chip_peak_file_dir = "../../../../data/wang_et_al/processed/run_1_ensembl52/"
atac_peak_file_dir = "../../../../data/wang_et_al/processed/run_2_ensembl52/"
sig_region_file = "./data/sig_regions.csv"

# Genome data
# gff_source = "../../../../data/AX4/genome_ver_2_7/gff3"
# chrom_lengths_file = "../../../../data/AX4/genome_ver_2_7/chromosome_lengths.txt"
# dictybase_cds_id_file = "../../../../data/AX4/genome_ver_2_7/fastas/dicty_primary_cds_ids.txt"
gff_source = "../../../../data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../../../data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
paralog_file = "./data/filtered/paralog_filt.tsv"
ensembl_cds_id_file = "../../../../data/AX4/genome_ver_2_7/ensembl_52/cds_ids.txt"
blacklist_file = "./data/blacklist_with_tes.csv"
final_gene_list = "./data/filtered/final_gene_list.txt"
singleton_file = "./data/filtered/singleton_filt.tsv"

# Save directory
ser_data_dir = "./data/julia_serialized/"

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
# peak_files = [filter(fn -> contains(fn, "FA"), peak_files); filter(fn -> contains(fn, "MA"), peak_files); filter(fn -> contains(fn, "VA"), peak_files)]
peak_data = binpeaks(peak_files, chrom_lengths_file)
addtogenes!(ref_genome, peak_data)

# Load singleton data
singletons = CSV.read(singleton_file, DataFrame).GeneID

# Load paralog data
paralog_data = CSV.read(paralog_file, DataFrame)
select!(paralog_data, ["GeneID", "ParalogID", "dS"])
filter!(row -> row["dS"] <= 3, paralog_data)

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
# k27ac_inds = [1, 4, 7]
# k4me3_inds = [2, 5, 8]
# k9me3_inds = [3, 6, 9]
sample_inds_vec = [k27ac_inds, 
                   k4me3_inds, 
                   k9me3_inds, 
                   atac_inds]
# sample_inds_vec = [k27ac_inds, k4me3_inds, k9me3_inds]
global_means_vec = 
[mean([mean([mean(getsiginrange(gene, GeneRange(TSS(), TES(), -500, 500), sample_ind)) for sample_ind in sample_inds]) for gene in filtered_gene_list if gene.id in singletons]) for sample_inds in sample_inds_vec]

tss_enrich = plot_enrich_region(paralog_data, 
    filtered_paralog_list, 
    sample_inds_vec, 
    [GeneRange(TSS(), TSS(), -500, 0) for i in 1:4], 
    fold_change_over_mean=true, 
    global_means=global_means_vec,
    z_min=0,
    z_max=4,
    return_figs=true)

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
bar_plots, p_vals, means_vecs = plot_bar(paralog_data, filtered_paralog_list, sample_inds_vec, [GeneRange(TSS(), TES(), sig_region_df.Start[1], parse(Int, sig_region_df.End[1])), # K27ac
                                                                            GeneRange(TSS(), TES(), sig_region_df.Start[2], parse(Int, sig_region_df.End[2])), # K4me3
                                                                            GeneRange(TSS(), TSS(), sig_region_df.Start[3], 100), # K9me3
                                                                            GeneRange(TSS(), TES(), sig_region_df.Start[4], parse(Int, sig_region_df.End[4]))], # ATAC
                                                            global_means_vec, [0,4], true, true)
adj_p_vals = adjust(pvalue.(p_vals), Bonferroni())
serialize(joinpath(ser_data_dir, "bar_plots_dS.jls"), bar_plots)
serialize_to_json(joinpath(ser_data_dir, "means_vecs_ds.json"), means_vecs)

