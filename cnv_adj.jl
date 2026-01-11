include("prelude.jl")

using Serialization
using CategoricalArrays

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
expr_data_file = "../../dicty_data/filtered/expr_data_filt_kallisto_2_single.tsv"
expr_data = CSV.read(expr_data_file, DataFrame)
rename!(expr_data, ["GeneID", "V", "S", "M", "F"])
select!(expr_data, ["GeneID", "F", "M", "V"])

# Average gene expression counts across all life-cycle stages
insertcols!(expr_data, :Avg => mean.(eachrow(expr_data[:, 2:end])))
select!(expr_data, ["GeneID", "Avg"])

# Genome data

gff_source = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/chromosome_lengths.txt"
singleton_list = CSV.read("../../dicty_data/singleton_list_unfilt.txt", DataFrame)[:,1]

# Input count data
input_count_experiment = deserialize("../../dicty_data/julia_serialized/input_count_experiment.jls")
input_counts = only(input_count_experiment.samples)

# Load a reference genome
ref_genome = loadgenome(gff_file_dir, chrom_lengths_file)

# Calculate average input read count for singletons
singleton_avg = mean([avg_input_count(String(gene_id), ref_genome, input_counts) for gene_id in singleton_list if gene_id in ref_genome.genes[1]])

# Get the genes with expression data
adjusted_expr_ratios = [id in ref_genome.genes[1] ? max(avg_input_count(id, ref_genome, input_counts) / singleton_avg, 1) : missing for id in expr_data.GeneID]
expr_data.Avg = [ismissing(adjusted_expr_ratios[i]) ? missing : expr_data.Avg[i] / adjusted_expr_ratios[i] for i in eachindex(adjusted_expr_ratios)]
CSV.write("../../dicty_data/filtered/expr_data_filt_adj.tsv", expr_data, delim='\t')