include("prelude.jl")

using .RepeatUtils

# Files:
repeat_file = "../../dicty_data/rm_genomes/ensembl_52/full/results_v1/onecodetofindthemall/Dictyostelium_discoideum.dicty_2.7.dna.toplevel_aggregated_compat.csv"
gff_source = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
ensembl_cds_id_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/cds_ids.txt"
blacklist_file = "./blacklists/cds_blacklist_full.tsv"
te_id_file = "./blacklists/ensembl_te_ids_with_predicted.tsv"

# Load genome data
ref_genome = loadgenome(gff_source, chrom_lengths_file)
numbered_chroms = ["$i" for i in 1:6]


# Add repeats to the reference genome object:
repeat_df = CSV.read(repeat_file, DataFrame)
convert_to_repeats!(ref_genome, repeat_df, allow_missing_scaffolds=true)
# addrepeats!(ref_genome, repeat_file)

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

# Caclulate the distance to the nearest repeat for each CDS gene. If no repeat is found on the same scaffold, 
# the distance is set to Inf.

repeat_dist_df = DataFrame(GeneID = [gene.id for gene in cds_genes],
    Scaffold = [gene.scaffold.name for gene in cds_genes],
    Distance = Vector{Union{Missing, Float64}}(zeros(Missing, length(cds_genes))),
    Repeat_start = Vector{Union{Missing, Float64}}(zeros(Missing, length(cds_genes))),
    Repeat_end = Vector{Union{Missing, Float64}}(zeros(Missing, length(cds_genes))))


for (i, gene) in enumerate(cds_genes)

    temp_scaffold = gene.scaffold
    repeat_dist_df[i, :GeneID] = gene.id
    repeat_dist_df[i, :Scaffold] = temp_scaffold.name
    distance = Inf
    repeat_start = Inf
    repeat_end = Inf

    for repeat_elem in temp_scaffold.repeats

        repeat_start = repeat_elem.repeat_start
        repeat_end = repeat_elem.repeat_end

        if GenomeTypes.hasoverlap(repeat_elem.repeat_start, gene.gene_start, repeat_elem.repeat_end, gene.gene_end)

            distance = 0
            break
        else

            distance = min(distance, min(abs(repeat_elem.repeat_start - gene.gene_end), abs(gene.gene_start - repeat_elem.repeat_end)))
        end
    end

    repeat_dist_df[i, :Distance] = distance
    repeat_dist_df[i, :Repeat_start] = repeat_start
    repeat_dist_df[i, :Repeat_end] = repeat_end
end

CSV.write("../../dicty_data/repeat_distance.csv", repeat_dist_df)