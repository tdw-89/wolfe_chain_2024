# Run third

using Serialization
using CSV
using DataFrames
using StatsBase
using CategoricalArrays

# Custom lib src:
include("./custom_lib/genomic_data.jl")
include("./custom_lib/load_gff.jl")
include("./custom_lib/misc_utils.jl")

function avg_input_count(gene_id::S, ref_genome::RefGenome, input_counts::SampleData) where S <: AbstractString
    gene = get(ref_genome, gene_id)
    gene_chrom = gene.scaffold.name
    gene_start = gene.gene_start
    gene_end = gene.gene_end
    chrom_ind = findfirst(chrom -> chrom.name == gene_chrom, input_counts.chroms)

    # return isnothing(chrom_ind) # DEBUG REMOVE

    return mean(input_counts.chroms[chrom_ind].signal[gene_start:gene_end])
end

# Expression data
expr_data_file = "./data/filtered/expr_data_filt_kallisto_2_single.tsv"
expr_data = CSV.read(expr_data_file, DataFrame)
rename!(expr_data, ["GeneID", "V", "S", "M", "F"])
select!(expr_data, ["GeneID", "F", "M", "V"])

# Average gene expression counts across all life-cycle stages
insertcols!(expr_data, :Avg => mean.(eachrow(expr_data[:, 2:end])))
select!(expr_data, ["GeneID", "Avg"])

# Genome data
gff_file_dir = "../../../../data/AX4/genome_ver_2_7/gff3"
chrom_lengths_file = "../../../../data/AX4/genome_ver_2_7/chromosome_lengths.txt"
singleton_list = CSV.read("./data/singleton_list_unfilt.txt", DataFrame)[:,1]

# Input count data
input_count_experiment = deserialize("./data/julia_serialized/input_count_experiment.jls")
input_counts = only(input_count_experiment.samples)

# Load a reference genome
ref_genome = loadgenome(gff_file_dir, chrom_lengths_file)

# Calculate average input read count for singletons
singleton_avg = mean([avg_input_count(String(gene_id), ref_genome, input_counts) for gene_id in singleton_list if gene_id in ref_genome.genes[1]])

# Get the genes with expression data
adjusted_expr_ratios = [id in ref_genome.genes[1] ? max(avg_input_count(id, ref_genome, input_counts) / singleton_avg, 1) : missing for id in expr_data.GeneID]
expr_data.Avg = [ismissing(adjusted_expr_ratios[i]) ? missing : expr_data.Avg[i] / adjusted_expr_ratios[i] for i in eachindex(adjusted_expr_ratios)]
CSV.write("./data/filtered/expr_data_filt_adj.tsv", expr_data, delim='\t')

# Old method:
#=
paralog_file = "./data/filtered/paralog_filt.tsv"

# Input count data
input_count_experiment = deserialize("./data/julia_serialized/input_count_experiment.jls")
input_counts = only(input_count_experiment.samples)

# Load a reference genome
ref_genome = loadgenome(gff_file_dir, chrom_lengths_file)

# Load paralog data
paralog_df = CSV.read(paralog_file, DataFrame)
paralog_ids = vcat(paralog_df[:,1], paralog_df[:,2])
paralog_ds_quants = levelcode.(cut(paralog_df.dS, 10))
paralod_df_low_ds = paralog_df[paralog_ds_quants .<= 2, :]
paralog_ids_low = vcat(paralod_df_low_ds[:,1], paralod_df_low_ds[:,2])

singleton_avg = log2(mean([avg_input_count(String(gene_id), ref_genome, input_counts) for gene_id in singleton_list if gene_id in ref_genome.genes[1]]))
singleton_ratio = [(gene_id, max(log2(avg_input_count(String(gene_id), ref_genome, input_counts)) / singleton_avg, 1)) for gene_id in singleton_list if gene_id in ref_genome.genes[1]]
singleton_ratio_df = DataFrame(singleton_ratio, [:GeneID, :InputRatio])

paralog_avg = mean([avg_input_count(String(gene_id), ref_genome, input_counts) for gene_id in paralog_ids])
paralog_avg_low = mean([avg_input_count(String(gene_id), ref_genome, input_counts) for gene_id in paralog_ids_low])

paralog_ratio_1 = [max(log2(avg_input_count(String(gene_id), ref_genome, input_counts)) / singleton_avg, 1) for gene_id in paralog_df[:,1]]
paralog_ratio_2 = [max(log2(avg_input_count(String(gene_id), ref_genome, input_counts)) / singleton_avg, 1) for gene_id in paralog_df[:,2]]
paralog_avg_ratio = mean(hcat(paralog_ratio_1, paralog_ratio_2), dims=2)[:,1]
paralog_avg_ratio_df = hcat(paralog_df, DataFrame(:InputRatio => paralog_avg_ratio))

CSV.write("./data/filtered/paralog_input_ratio.tsv", paralog_avg_ratio_df, delim='\t')
CSV.write("./data/filtered/singleton_input_ratio.tsv", singleton_ratio_df, delim='\t')
=#