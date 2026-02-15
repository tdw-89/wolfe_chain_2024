include("prelude.jl")

using .RepeatUtils

# Files:
te_id_file = "blacklists/ensembl_te_ids_with_predicted.tsv"
pred_te_coord_file = "../../dicty_data/rm_genomes/ensembl_52/full/results_v1/onecodetofindthemall/Dictyostelium_discoideum.dicty_2.7.dna.toplevel_aggregated_compat.csv"
gff_source = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
chip_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_1_ensembl52/"

peak_files = [filter(fn -> !contains(fn, r"_S[AB]_"), readdir(chip_peak_file_dir, join=true)); readdir(atac_peak_file_dir, join=true)]

# Load genome data
ref_genome = loadgenome(gff_source, chrom_lengths_file)

te_ids = open(te_id_file) do file
    readlines(file)
end

# Load the predicted TE coordinates:
pred_te_df = CSV.read(pred_te_coord_file, DataFrame)

# Add the predicted TEs to the reference genome object:
convert_to_repeats!(ref_genome, pred_te_df, allow_missing_scaffolds=true)

# Convert TE genes to 'repeats' in the reference genome object:
te_genes_list = [gene for gene in ref_genome.genes[2] if gene.id in te_ids]
move_to_repeats!(ref_genome, te_genes_list)
te_genes_list = nothing