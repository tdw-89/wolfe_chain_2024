include("prelude.jl")

using CategoricalArrays
using Serialization

LOW_DS_QUANT_NUM_CUTOFF = 4

using .EnrichmentUtils
using .MiscUtils

human_gff = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.99.gff3"
human_paralog_info = "../../dicty_data/filtered/human_paralog_info_filt.csv"
chrom_lengths_file = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/chromosome_lengths.txt"
singleton_list_file = "../../dicty_data/filtered/human_singletons_filt.csv"
paralog_file = "../../dicty_data/filtered/human_paralog_info_filt.csv"

# Load the reference
ref_genome = loadgenome(human_gff, chrom_lengths_file)

# Load the peak data
peak_data = deserialize("../../dicty_data/julia_serialized/human_h3k9me3_exper.jls")

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

quant_ids = vcat([String.(paralog_data.GeneID[quantile_vals .== q]) for q in unique(sort(quantile_vals))])

low_ds_ids = String.(vcat(paralog_data.GeneID[quantile_vals .<= LOW_DS_QUANT_NUM_CUTOFF],
                          paralog_data.ParalogID[quantile_vals .<= LOW_DS_QUANT_NUM_CUTOFF]))
high_ds_ids = String.(vcat(paralog_data.GeneID[quantile_vals .> LOW_DS_QUANT_NUM_CUTOFF],
                           paralog_data.ParalogID[quantile_vals .> LOW_DS_QUANT_NUM_CUTOFF]))
low_ds_paralog_genes = get(ref_genome, low_ds_ids)
high_ds_paralog_genes = get(ref_genome, high_ds_ids)

low_ds_paralog_has_K9me3 = [
    sum([
        sum(getsiginrange(gene, GeneRange(TSS(), TES(), -500, 500), ind)) > 1
     for ind in eachindex(gene.samples)
    ]) > 1
    for gene in low_ds_paralog_genes
]

high_ds_paralog_has_K9me3 = [
    sum([
        sum(getsiginrange(gene, GeneRange(TSS(), TES(), -500, 500), ind)) > 1
     for ind in eachindex(gene.samples)
    ]) > 1
    for gene in high_ds_paralog_genes
]

low_ds_n_has_k9me3 = sum(low_ds_paralog_has_K9me3)
low_ds_n_no_k9me3 = length(low_ds_paralog_has_K9me3) - low_ds_n_has_k9me3
high_ds_n_has_K9me3 = sum(high_ds_paralog_has_K9me3)
high_ds_n_no_K9me3 = length(high_ds_paralog_has_K9me3) - high_ds_n_has_K9me3

cont_table_low_vs_other = [low_ds_n_has_k9me3 high_ds_n_has_K9me3;
                           low_ds_n_no_k9me3 high_ds_n_no_K9me3]

pvalue(FisherExactTest(cont_table_low_vs_other[1, 1], cont_table_low_vs_other[1, 2],
                       cont_table_low_vs_other[2, 1], cont_table_low_vs_other[2, 2]))