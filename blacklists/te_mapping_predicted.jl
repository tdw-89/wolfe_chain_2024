#=
STEP 2: Get a list of TE IDs from the DB reference and those overlapping a predicted TE
=#

include("../prelude.jl")
# Given that the numbered chromosomes were the same in the dictybase release and Ensembly releases,
# only the non-numbered (contigs) chromosomes were re-scanned for TEs with RepeatModeler

using CSV
using DataFrames

using .RepeatUtils

cd(@__DIR__)

# Genome data:
gff_source = "../../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
paralog_file = "../../../dicty_data/filtered/paralog_filt.tsv"

# Predicted TE files (Ensembl 52):
tes_ensembl = "../../../dicty_data/rm_genomes/ensembl_52/full/results_v1/onecodetofindthemall/Dictyostelium_discoideum.dicty_2.7.dna.toplevel_aggregated_compat.csv"

# mapped TEs file:
mapped_tes_file = "ensembl_te_ids_full.tsv"

# Load the genome data:
ref_genome = loadgenome(gff_source, chrom_lengths_file)
paralog_data = CSV.read(paralog_file, DataFrame)

# Load the TE data:
te_df_ensembl = CSV.read(tes_ensembl, DataFrame)
mapped_tes = CSV.read(mapped_tes_file, DataFrame)

# Add the repeats to the genome:
convert_to_repeats!(ref_genome, te_df_ensembl; allow_missing_scaffolds=true)

# Get genes that overlap with predicted TEs:
te_genes = [findoverlappinggenes(repeat_elem) for repeat_elem in ref_genome.repeats]
overlapping_ids = [te_genes[i][1] for i in eachindex(te_genes) if !ismissing(te_genes[i])]
overlapping_ids = reduce(vcat, overlapping_ids)
overlapping_ids = unique(overlapping_ids)

# Merge the overlapping IDs with the TE IDs lifted from the dictybase genome:
all_te_IDs = union(mapped_tes.GeneID, overlapping_ids)
CSV.write("ensembl_te_ids_with_predicted.tsv", DataFrame(GeneID=all_te_IDs))