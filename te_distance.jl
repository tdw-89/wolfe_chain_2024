include("prelude.jl")

using FastaIO

using .RepeatUtils

# Files:
te_id_file = "blacklists/ensembl_te_ids_with_predicted.tsv"
pred_te_coord_file = "../../dicty_data/rm_genomes/ensembl_52/full/results_v1/onecodetofindthemall/Dictyostelium_discoideum.dicty_2.7.dna.toplevel_aggregated_compat.csv"
gff_source = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
ensembl_cds_id_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/cds_ids.txt"
blacklist_file = "./blacklists/cds_blacklist_full.tsv"

# Load genome data
ref_genome = loadgenome(gff_source, chrom_lengths_file)

# Get CDS genes not in blacklist and not in TE list:
cds_ids = open(ensembl_cds_id_file) do file
    readlines(file)
end

te_ids = open(te_id_file) do file
    readlines(file)
end

blacklist_ids = open(blacklist_file) do file
    readlines(file)
end

blacklist_ids = blacklist_ids[2:end]

cds_genes = [gene for gene in ref_genome.genes[2] if gene.id in cds_ids && gene.id ∉ te_ids && gene.id ∉ blacklist_ids]

# Load the predicted TE coordinates:
pred_te_df = CSV.read(pred_te_coord_file, DataFrame)

# Add the predicted TEs to the reference genome object:
convert_to_repeats!(ref_genome, pred_te_df, allow_missing_scaffolds=true)

# Convert TE genes to 'repeats' in the reference genome object:
te_genes_list = [gene for gene in ref_genome.genes[2] if gene.id in te_ids]
move_to_repeats!(ref_genome, te_genes_list)
te_genes_list = nothing

# Caclulate the distance to the nearest TE for each CDS gene. If no TE is found on the same scaffold, 
# the distance is set to Inf.
te_dist_df = DataFrame(GeneID = [gene.id for gene in cds_genes],
    Scaffold = [gene.scaffold.name for gene in cds_genes],
    Distance = Vector{Union{Missing, Float64}}(zeros(Missing, length(cds_genes))),
    TE_start = Vector{Union{Missing, Float64}}(zeros(Missing, length(cds_genes))),
    TE_end = Vector{Union{Missing, Float64}}(zeros(Missing, length(cds_genes))))

for (i, gene) in enumerate(cds_genes)

    temp_scaffold = gene.scaffold
    te_dist_df[i, :GeneID] = gene.id
    te_dist_df[i, :Scaffold] = temp_scaffold.name
    distance = Inf
    te_start = Inf
    te_end = Inf

    for repeat_elem in temp_scaffold.repeats

        te_start = repeat_elem.repeat_start
        te_end = repeat_elem.repeat_end

        if GenomeTypes.hasoverlap(repeat_elem.repeat_start, gene.gene_start, repeat_elem.repeat_end, gene.gene_end)

            distance = 0
            break
        else

            distance = min(distance, 
                       min(abs(repeat_elem.repeat_start - gene.gene_end), 
                           abs(gene.gene_start - repeat_elem.repeat_end)))
        end
    end

    te_dist_df[i, :Distance] = distance
    te_dist_df[i, :TE_start] = te_start
    te_dist_df[i, :TE_end] = te_end
end

CSV.write("../../dicty_data/te_distance.csv", te_dist_df)