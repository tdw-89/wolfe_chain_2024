using CSV
using DataFrames
using FastaIO

# custom lib:
include("custom_lib/load_gff.jl")
include("custom_lib/te_utils.jl")
include("custom_lib/misc_utils.jl")

# Functions:

# function parse_rmsk_info(rmsk_info)
#     """
#     Parse the rmsk info field.
#     """
#     rmsk_info = strip.(split(rmsk_info, ";"))
#     family_info = rmsk_info[3] # NOTE: All strings from column 9 were checked and have family as the 3rd element and class as the 4th in the transcript TE .gtf file from Hammell lab
#     class_info = rmsk_info[4]
#     return split(family_info, r"\"")[2], 
#             split(class_info, r"\"")[2]
# end

# Files
human_gff = "../../../../data/mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.99.gff3"
chrom_lengths_file = "../../../../data/mammals/primates/h_sapiens/Ensembl_99/chromosome_lengths.txt"
repeat_file = "../../../../data/mammals/primates/h_sapiens/GRCh38_UCSC_rmsk.tsv.gz" # UCSC-derived table of repeatmasker-identified repeats
ensembl_cds_file = "../../../../data/mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.cds.all.fa.gz"

# Load the reference
human_ref = loadgenome(human_gff, chrom_lengths_file)

# Load the repeat df
repeat_df = CSV.read(repeat_file, DataFrame)

# Filter out rna genes
filter!(row -> !contains(row.repClass, "RNA"), repeat_df)

# Select and rename
select!(repeat_df, [6, 7, 8, 12, 13])
rename!(repeat_df, [:Chromosome, :Start, :End, :Type, :Family])
repeat_df.Chromosome = map(chr -> String(split(chr, "chr")[2]), repeat_df.Chromosome)

# Add repeats to the reference genome
convert_to_repeats!(human_ref, repeat_df, allow_missing_scaffolds=true)

# Get the CDS IDs
cds_ids = String[]
FastaReader(ensembl_cds_file) do fr
    for (name, seq) in fr
        if occursin("gene_biotype:protein_coding", name) # Excludes T-cell receptor genes, Imunnoglobin genes, and 'processed_pseudogenes'
            push!(cds_ids, match(r"gene:[A-Za-z0-9]+", name).match[6:end])
        
        end
    end
end

unique!(cds_ids)
@assert all(length.(cds_ids) .== 15) && all(contains.(cds_ids, "ENSG"))

# Collect the CDS genes from the reference object
cds_genes = get(human_ref, cds_ids)

# Identify the genes from the CDS file missing from the reference object
missing_genes = cds_ids[findall(ismissing, cds_genes)] # All missing genes have chromosome info with the following string: "chromosome:GRCh38:CHR_H", which are
                                                       # all patch chromosomes/sequences

# Collect the non-missing genes
cds_genes = collect(skipmissing(cds_genes))

# Create the repeat distance df
repeat_dist_df = DataFrame(GeneID = [gene.id for gene in cds_genes],
    Scaffold = [gene.scaffold.name for gene in cds_genes],
    Distance = Vector{Union{Missing, Float64}}(zeros(Missing, length(cds_genes))),
    TE_start = Vector{Union{Missing, Float64}}(zeros(Missing, length(cds_genes))),
    TE_end = Vector{Union{Missing, Float64}}(zeros(Missing, length(cds_genes))))


progress_bits = BitVector(zeros(Bool, length(cds_genes)))
perc_increment = length(cds_genes) ÷ 100
Threads.@threads for i in eachindex(cds_genes)

    temp_gene = cds_genes[i]
    temp_scaffold = temp_gene.scaffold
    repeat_dist_df[i, :GeneID] = temp_gene.id
    repeat_dist_df[i, :Scaffold] = temp_scaffold.name
    distance = Inf
    repeat_start = Inf
    repeat_end = Inf

    for repeat_elem in temp_scaffold.repeats

        repeat_start = repeat_elem.repeat_start
        repeat_end = repeat_elem.repeat_end

        if hasoverlap(repeat_elem.repeat_start, temp_gene.gene_start, repeat_elem.repeat_end, temp_gene.gene_end)

            distance = 0
            break
        else

            distance = min(distance, min(abs(repeat_elem.repeat_start - temp_gene.gene_end), abs(temp_gene.gene_start - repeat_elem.repeat_end)))
        end
    end

    repeat_dist_df[i, :Distance] = distance
    repeat_dist_df[i, :TE_start] = repeat_start
    repeat_dist_df[i, :TE_end] = repeat_end
    progress_bits[i] = 1

    if sum(progress_bits) % perc_increment == 0
        println("Progress: ", sum(progress_bits) ÷ perc_increment, "%")

    end
end

# Save the distance df
CSV.write("./data/repeat_distance_human.csv", repeat_dist_df)
CSV.write("./data/repeat_distance_human_numbered_chroms.csv", filter(row -> row.Scaffold in [["$i" for i in 1:22]; "X"; "Y"], repeat_dist_df))