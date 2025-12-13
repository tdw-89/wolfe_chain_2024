using CSV
using DataFrames
using Serialization
using CategoricalArrays
using HypothesisTests
using JSON

include("custom_lib/load_gff.jl")
include("custom_lib/genomic_data.jl")
include("custom_lib/enrichment_utils.jl")
include("custom_lib/te_utils.jl")

reload_peak_data = false

# function for serializing enrichment vectors to JSON for use in R
function serialize_to_json(file_path, vecs)
    mark_names = ["K9me3"]
    
    open(file_path, "w") do file
        JSON.print(file, Dict([mark => vs for (mark, vs) in zip(mark_names, vecs)]), 4)

    end
end

# Human genomic data
human_gff = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.99.gff3"
human_paralog_info = "../../dicty_data/filtered/human_paralog_info_filt.csv"
peak_data_dir = "../../dicty_data/mammals/primates/h_sapiens/ENCODE_histone_mods/"
chrom_lengths_file = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/chromosome_lengths.txt"
cds_id_file = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/CDS_IDs.txt"
nmd_id_file = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/nmd_candidate_ids.txt"
singleton_file = "../../dicty_data/filtered/human_singletons_filt.csv"

# Load the reference
human_ref = loadgenome(human_gff, chrom_lengths_file)
cds_df = CSV.read(cds_id_file, DataFrame, header=false)
nmd_df = CSV.read(nmd_id_file, DataFrame, header=false)

# Add the peak data
if reload_peak_data
    peak_files = reduce(vcat, [[joinpath(root, fn) for fn in files] for (root, dirs, files) in walkdir(peak_data_dir)])
    filter!(fn -> contains(fn, ".bed"), peak_files)
    peak_data = binpeaks(peak_files, chrom_lengths_file)
    serialize("../../dicty_data/julia_serialized/human_h3k9me3_exper.jls", peak_data)

else
    peak_data = deserialize("../../dicty_data/julia_serialized/human_h3k9me3_exper.jls")

end

addtogenes!(human_ref, peak_data)
peak_data = nothing
GC.gc()

# Load the paralog info
paralog_data = CSV.read(human_paralog_info, DataFrame)
filter!(row -> row.dS <= 3, paralog_data)
CSV.write("data/paralog_data_human.csv", paralog_data)
select!(paralog_data, ["GeneID", "ParalogID", "dS"])


# Get the global mean enrichment for H3K9me3 from TSS-500:TSS+500 for all coding genes
global_mean = mean([
                mean([mean(getsiginrange(gene, GeneRange(TSS(), TES(), -500, 500), i))
                for i in eachindex(human_ref.genes[2][1].samples)])
                for gene in human_ref.genes[2]
                if gene.id in cds_df[!,1] && gene.id ∉ nmd_df[!,1]
])

# Plot the enrichment
tss_enrich = plot_enrich_region(paralog_data,
                                human_ref.genes[2],
                                [collect(eachindex(human_ref.genes[2][1].samples))],
                                [GeneRange(TSS(), TSS(), -500, -1)],
                                fold_change_over_mean=true,
                                global_means=[global_mean],
                                z_min=0,
                                z_max=4,
                                return_figs=true)

serialize("../../dicty_data/julia_serialized/human_tss_enrich_plots_dS.jls", tss_enrich)

body_enrich = plot_enrich_percent(paralog_data,
                    human_ref.genes[2],
                    [collect(eachindex(human_ref.genes[2][1].samples))],
                    fold_change_over_mean=true,
                    global_means=[global_mean],
                    z_min=0,
                    z_max=4,
                    return_figs=true)
serialize("../../dicty_data/julia_serialized/human_body_enrich_plots_dS.jls", body_enrich)

tes_enrich = plot_enrich_region(paralog_data,
                                human_ref.genes[2],
                                [collect(eachindex(human_ref.genes[2][1].samples))],
                                [GeneRange(TES(), TES(), 1, 500)],
                                fold_change_over_mean=true,
                                global_means=[global_mean],
                                z_min=0,
                                z_max=4,
                                return_figs=true)
serialize("../../dicty_data/julia_serialized/human_tes_enrich_plots_dS.jls", tes_enrich)

bar_plots, p_vals, means_vecs = plot_bar(paralog_data,
                    human_ref.genes[2],
                    [collect(eachindex(human_ref.genes[2][1].samples))],
                    [GeneRange(TSS(), TES(), -500, 500)],
                    [global_mean],
                    [0,4],
                    true,
                    true);

p_vals_perm_cor = get_cor(paralog_data,
                          GeneRange(TSS(),TES(), -500, 500),
                          [1:16...],
                          human_ref,
                          global_mean)

serialize("../../dicty_data/julia_serialized/human_bar_plots_dS.jls", bar_plots)
serialize_to_json("../../dicty_data/julia_serialized/means_vecs_dS_human.json", means_vecs)