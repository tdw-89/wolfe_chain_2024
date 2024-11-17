using HypothesisTests
using CSV
using DataFrames
using CategoricalArrays
using PlotlyJS
using Serialization

include("./custom_lib/load_gff.jl")
include("./custom_lib/enrichment_utils.jl")
include("./custom_lib/misc_utils.jl")

human_gff = "../../../../data/mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.99.gff3"
human_paralog_info = "./data/filtered/human_paralog_info_filt.csv"
chrom_lengths_file = "../../../../data/mammals/primates/h_sapiens/Ensembl_99/chromosome_lengths.txt"
singleton_list_file = "./data/filtered/human_singletons_filt.csv"
paralog_file = "./data/filtered/human_paralog_info_filt.csv"

# Load the reference
ref_genome = loadgenome(human_gff, chrom_lengths_file)

# Load the peak data
peak_data = deserialize("./data/julia_serialized/human_h3k9me3_exper.jls")

# Add the peak data to the reference:
addtogenes!(ref_genome, peak_data)

# Load paralog and singleton data
paralog_data = CSV.read(paralog_file, DataFrame)
singleton_list = CSV.read(singleton_list_file, DataFrame)

# Filter paralogs according to id threshold
select!(paralog_data, ["GeneID", "ParalogID", "dS"])
filter!(row -> row.dS <= 3.0, paralog_data)

# Get the dS values and quantiles for the filtered pairs
quantile_labels = cut(paralog_data[!, 3], 10)
quantile_vals = levelcode.(quantile_labels)

low_ds_ids = String.(vcat(paralog_data.GeneID[quantile_vals .<= 4], 
                          paralog_data.ParalogID[quantile_vals .<= 4]))
paralog_genes = get(ref_genome, low_ds_ids)

paralog_gene_has_K9me3 = [
    siginrange(gene, GeneRange(TSS(), TES(), -500, 500), 1) ?
    sum([
        sum(getsiginrange(gene, GeneRange(TSS(), TES(), -500, 500), ind)) > 1
     for ind in eachindex(gene.samples)
    ]) > 1 :
    sum([
        sum(getsiginrange(gene, GeneRange(TSS(), TES(), 0, 0), ind)) > 1
     for ind in eachindex(gene.samples)
    ]) > 1
    for gene in paralog_genes
]

low_ds_n_has_k9me3 = sum(paralog_gene_has_K9me3)
n_total = length(paralog_genes)

singleton_genes = [gene for gene in ref_genome.genes[2] if gene.id in singleton_list.GeneID]