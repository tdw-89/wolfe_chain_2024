include("prelude.jl")

using .EnrichmentUtils

using Serialization
using JSON
using Random

TISSUE_TYPE = "All"  # Options: "All", "V", "M", "F"

function get_sample_inds(tissue_type)
    if tissue_type == "All"
        return [
            [1, 4, 7], # H3K27ac
            [2, 5, 8], # H3K4me3
            [3, 6, 9], # H3K9me3
            [10, 11, 12] # ATAC
            ]
    elseif tissue_type == "V"
        return [[1], [2], [3]]
    elseif tissue_type == "M"
        return [[1], [2], [3], [4]]
    elseif tissue_type == "F"
        return [[1], [2], [3], [4]]
    else
        error("Invalid TISSUE_TYPE: $tissue_type")
    end
end

# function for serializing enrichment vectors to JSON for use in R
function serialize_to_json(file_path, vecs)
    mark_names = ["K27ac", "K4me3", "K9me3", "ATAC"]
    
    open(file_path, "w") do file
        JSON.print(file, Dict([mark => vs for (mark, vs) in zip(mark_names, vecs)]), 4)

    end
end

gene_range =
    GeneRange(
    TSS(),
    TES(),
    -500,
    500
    )

# Peak files
chip_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_1_ensembl52/"
atac_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_2_ensembl52/"

# Genome data
gff_source = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
paralog_file = "../../dicty_data/filtered/paralog_filt.tsv"
ensembl_cds_id_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/cds_ids.txt"
blacklist_file = "../../dicty_data/blacklist_with_tes.csv"
final_gene_list = "../../dicty_data/filtered/final_gene_list.txt"

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

# NOTE: If 'V' then there will be no ATAC data, if 'S'
if TISSUE_TYPE != "All"
    peak_files = filter(fn -> contains(fn, "_$TISSUE_TYPE"), peak_files)
end

peak_data = binpeaks(peak_files, chrom_lengths_file)
addtogenes!(ref_genome, peak_data)

# Load paralog data
paralog_data = CSV.read(paralog_file, DataFrame)
filter!(row -> row["dS"] <= 3, paralog_data)
CSV.write("../../dicty_data/filtered/paralog_ds_filt.csv", paralog_data)
fig = plot(histogram(x=paralog_data.dS), Layout(title="Distribution of dS values for paralog pairs", xaxis_title="dS", yaxis_title="Count"))
savefig(fig, "../../dicty_data/filtered/paralog_ds_hist.html")
select!(paralog_data, ["GeneID", "ParalogID", "dS"])

# paralog_data[findall(map(id -> !siginrange(get(ref_genome, id),  GeneRange(TSS(), TES(), -500, 500)), paralog_data.ParalogID)),:]

# select!(paralog_data, ["GeneID", "ParalogID", "MaxPerc"])
paralog_ids = vcat(paralog_data[:, 1], paralog_data[:, 2])

# Filter to genes that have 500 bp upstream and 250 downstream before a chromosome end
filtered_gene_list = [gene for gene in ref_genome.genes[2] if gene.id in final_gene_list]
filtered_paralog_list = [gene for gene in filtered_gene_list if gene.id in paralog_ids]

# Plot enrichment vs expression
sample_inds_vec = get_sample_inds(TISSUE_TYPE)
# get the global mean enrichment for each sample across all filtered genes
get_mean_sig(gene, gene_range, sample_inds) = mean([mean(getsiginrange(gene, gene_range, sample_ind)) for sample_ind in sample_inds])
get_global_mean(gene_list, sample_inds) = mean([get_mean_sig(gene, gene_range, sample_inds) for gene in gene_list])
get_global_std_dev(gene_list, sample_inds) = std([get_mean_sig(gene, gene_range, sample_inds) for gene in gene_list])
get_sig_dist(gene_list, sample_inds) = 
    EnrichmentUtils.SignalDistribution(
        get_global_mean(gene_list, sample_inds), 
        get_global_std_dev(gene_list, sample_inds)
    )
global_means_vec = [get_sig_dist(filtered_gene_list, sample_inds) for sample_inds in sample_inds_vec]

tss_enrich = plot_enrich_region(
    paralog_data, 
    filtered_paralog_list,
    sample_inds_vec,
    [GeneRange(TSS(), TSS(), -500, -1) for i in 1:length(sample_inds_vec)],
    fold_change_over_mean=false,
    global_means=global_means_vec,
    z_min=z_min,
    z_max=z_max,
    return_figs=true
    )

serialize(joinpath(ser_data_dir, "tss_enrich_plots_dS_$TISSUE_TYPE.jls"), tss_enrich)

body_enrich = plot_enrich_percent(paralog_data,
filtered_paralog_list, 
    sample_inds_vec, 
    fold_change_over_mean=false, 
    global_means=global_means_vec,
    z_min=z_min,
    z_max=z_max,
    return_figs=true)

serialize(joinpath(ser_data_dir, "body_enrich_plots_dS_$TISSUE_TYPE.jls"), body_enrich)

tes_enrich = plot_enrich_region(paralog_data,
filtered_paralog_list,
    sample_inds_vec,
    [GeneRange(TES(), TES(), 1, 500) for i in 1:length(sample_inds_vec)],
    fold_change_over_mean=false,
    global_means=global_means_vec,
    z_min=z_min,
    z_max=z_max,
    return_figs=true)

serialize(joinpath(ser_data_dir, "tes_enrich_plots_dS_$TISSUE_TYPE.jls"), tes_enrich)

# Plot bar plots of coverage in significant regions
bar_plots, kw_tests, means_vecs = 
    plot_bar(
        paralog_data,
        filtered_paralog_list, 
        sample_inds_vec,
        [gene_range for i in 1:length(sample_inds_vec)],
        global_means_vec, 
        [z_min, z_max], 
        fold_change_over_mean=false, 
        return_figs=true,
        horizontal=false,
        box_plots=true);
eta_sq_vals = η².(kw_tests)

p_vals_perm_cor = [get_cor(paralog_data, gene_range, sample_ind, ref_genome, global_mean.mean) for (sample_ind, gene_range, global_mean) in zip(sample_inds_vec,
                            [gene_range for i in 1:length(sample_inds_vec)],
                             global_means_vec)]
adj_p_vals_perm_cor = adjust([pair[1] for pair in p_vals_perm_cor], BenjaminiHochberg())
adj_p_vals_kw = adjust(pvalue.(kw_tests), BenjaminiHochberg())
serialize(joinpath(ser_data_dir, "bar_plots_dS_$TISSUE_TYPE.jls"), bar_plots)
serialize_to_json(joinpath(ser_data_dir, "means_vecs_ds_$TISSUE_TYPE.json"), means_vecs)