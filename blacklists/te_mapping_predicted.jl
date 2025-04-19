using CSV
using DataFrames

cd("blacklists")

# Custom lib src:
include("../custom_lib/load_gff.jl")
include("../custom_lib/te_utils.jl")

# Genome data:
gff_source = "../../../../../data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../../../../data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
chrom_names_file_db = "../../../../../data/AX4/genome_ver_2_7/chromosome_names.csv"
paralog_file = "../data/filtered/paralog_filt.tsv"

# Predicted TE files:

    # dictybase:
tes_dictybase = "../../../../../data/AX4/genome_ver_2_7/dicty_repeats-families.stk"

    # Ensembl 52:
tes_ensembl = "../../../../../data/AX4/genome_ver_2_7/repeat_modeler_output/RM_2530911.TueFeb132059322024/families-classified.stk"

# mapped TEs file:
mapped_tes_file = "ensembl_te_ids_full.tsv"

# Load the genome data:
chrom_names = CSV.read(chrom_names_file_db, DataFrame)
chrom_name_dict = Dict(chromid => chromname for (chromid, chromname) in zip(chrom_names[!, :ChromID], chrom_names[!, :ChromName]))
ref_genome = loadgenome(gff_source, chrom_lengths_file)
paralog_data = CSV.read(paralog_file, DataFrame)

# Load the repeat/TE data:
te_df_ensembl = parsestk(tes_ensembl)
te_df_dictybase = parsestk(tes_dictybase)
mapped_tes = CSV.read(mapped_tes_file, DataFrame)

# Filter to transposons only:
te_df_ensembl_filt = filter(row -> !ismissing(row.Type) && contains(row.Type, "Transposable"), te_df_ensembl)
te_df_dictybase_filt = filter(row -> !ismissing(row.Type) && contains(row.Type, "Transposable"), te_df_dictybase)

# Filter dictybase TEs to numbered chromosomes only:
te_df_dictybase_filt_chrom = filter(row -> chrom_name_dict[row.Chromosome] in ["$i" for i in 1:6], te_df_dictybase_filt)
te_df_dictybase_filt_chrom[!, :Chromosome] = [chrom_name_dict[row.Chromosome] for row in eachrow(te_df_dictybase_filt_chrom)]

# Combine the TE data:
te_df_combined = vcat(te_df_ensembl_filt, te_df_dictybase_filt_chrom)
CSV.write("predicted_te_combined.csv", te_df_combined)

# Add the repeats to the genome:
convert_to_repeats!(ref_genome, te_df_combined)

# Get genes that overlap with predicted TEs:
te_genes = [findoverlappinggenes(repeat_elem) for repeat_elem in ref_genome.repeats]
overlapping_ids = [te_genes[i][1] for i in eachindex(te_genes) if !ismissing(te_genes[i])]
overlapping_ids = reduce(vcat, overlapping_ids)
overlapping_ids = unique(overlapping_ids)

# Merge the overlapping IDs with the TE IDs lifted from the dictybase genome:
all_te_IDs = union(mapped_tes.GeneID, overlapping_ids)
CSV.write("ensembl_te_ids_with_predicted.tsv", DataFrame(GeneID=all_te_IDs))