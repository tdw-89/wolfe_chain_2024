include("prelude.jl")

using FastaIO
using Serialization

using .RepeatUtils
using .MiscUtils

only_h3k9me3 = false
only_specific_type = true
reload_peak_data = false
only_numbered_chroms = true
specific_type = "TE"

# Files
human_gff = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.99.gff3"
chrom_lengths_file = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/chromosome_lengths.txt"
repeat_file = "../../dicty_data/mammals/primates/h_sapiens/GRCh38_UCSC_rmsk.tsv.gz" # UCSC-derived table of repeatmasker-identified repeats
repeat_has_h3k9me3_file = "../../dicty_data/te_has_h3k9me3.csv"
ensembl_cds_file = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.cds.all.fa.gz"
peak_data_dir = "../../dicty_data/mammals/primates/h_sapiens/ENCODE_histone_mods/"

# Load the reference
human_ref = loadgenome(human_gff, chrom_lengths_file)
human_ref_2 = loadgenome(human_gff, chrom_lengths_file, feature_type="all");

# Load the repeat df
te_df = only_h3k9me3 ? CSV.read(repeat_has_h3k9me3_file, DataFrame) : CSV.read(repeat_file, DataFrame)

if !only_h3k9me3
    # Filter out low-confidence classes
    all_classes = filter(class -> !contains(class, "?"), unique(te_df.repClass))
    filter!(row -> row.repClass in all_classes, te_df)

    # Select and rename
    select!(te_df, [6, 7, 8, 12, 13])
    rename!(te_df, [:Chromosome, :Start, :End, :Type, :Family])
    te_df.Chromosome = map(chr -> String(split(chr, "chr")[2]), te_df.Chromosome)
    te_df.Start = te_df.Start .+ 1 # Convert to 1-based indexing

end

if only_specific_type
    if specific_type == "TE"
        te_types = ["LINE", 
                    "SINE", 
                    "DNA", 
                    "Retroposon", 
                    "LTR"]
        filter!(row -> row.Type in te_types, te_df)
    else
        if specific_type in te_df.Type
            filter!(row -> row.Type == specific_type, te_df)
        elseif specific_type in te_df.Family
            filter!(row -> row.Family == specific_type, te_df)
        end
    end
end

# Calculate total TE coverage statistics
if specific_type == "TE"
    dna_te_df = filter(row -> row.Type == "DNA", te_df)
    rte_te_df = filter(row -> row.Type != "DNA", te_df)
    dna_bases = sum(dna_te_df.End .- dna_te_df.Start .+ 1)
    rte_bases = sum(rte_te_df.End .- rte_te_df.Start .+ 1)
    total_te_bases = sum(te_df.End .- te_df.Start .+ 1)
    total_genome_bases = sum(values(human_ref.scaffolds) .|> s -> s.scaffold_end - s.scaffold_start + 1)
    
    dna_te_perc = dna_bases / total_genome_bases * 100.0
    rte_te_perc = rte_bases / total_genome_bases * 100.0
    total_te_perc = total_te_bases / total_genome_bases * 100.0
    println("DNA TE coverage: $dna_bases / $total_genome_bases * 100 = $dna_te_perc")
    println("Retrotransposon TE coverage: $rte_bases / $total_genome_bases * 100 = $rte_te_perc")
    println("Total TE coverage: $total_te_bases / $total_genome_bases * 100 = $total_te_perc")
end

# Filter to only include the main chromosomes
main_chroms = [["$i" for i in 1:22]..., "X", "Y"]
if only_numbered_chroms
    filter!(row -> row.Chromosome in main_chroms, te_df)

end

# Add repeats to the reference genome
convert_to_repeats!(human_ref, te_df, allow_missing_scaffolds=true)

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
missing_genes = cds_ids[findall(ismissing, cds_genes)] # All missing genes have chromosome info with the following string: "chromosome:GRCh38:CHR_H". This signifies that the gene is from an alternate haplotype

# Collect the non-missing genes
cds_genes = filter(gene -> gene.scaffold.name ∈ main_chroms, collect(skipmissing(cds_genes)))

# Create the repeat distance df
te_dist_df = DataFrame(GeneID = [gene.id for gene in cds_genes],
    Scaffold = [gene.scaffold.name for gene in cds_genes],
    Distance = Vector{Union{Missing, Float64}}(zeros(Missing, length(cds_genes))),
    TE_start = Vector{Union{Missing, Float64}}(zeros(Missing, length(cds_genes))),
    TE_end = Vector{Union{Missing, Float64}}(zeros(Missing, length(cds_genes))),
    Gene_start = [gene.gene_start for gene in cds_genes],
    Gene_end = [gene.gene_end for gene in cds_genes])

progress_bits = BitVector(zeros(Bool, length(cds_genes)))
perc_increment = length(cds_genes) ÷ 100

Threads.@threads for i in eachindex(cds_genes)
    temp_gene = cds_genes[i]
    temp_scaffold = temp_gene.scaffold
    te_dist_df[i, :GeneID] = temp_gene.id
    te_dist_df[i, :Scaffold] = temp_scaffold.name
    distance = Inf
    closest_repeat_start = Inf
    closest_repeat_end = Inf

    for repeat_elem in temp_scaffold.repeats
        repeat_start = repeat_elem.repeat_start
        repeat_end = repeat_elem.repeat_end

        if GenomeTypes.hasoverlap(repeat_elem.repeat_start, 
                        temp_gene.gene_start, 
                        repeat_elem.repeat_end, 
                        temp_gene.gene_end)
            distance = 0
            break

        else
            temp_distance = min(abs(repeat_elem.repeat_start - temp_gene.gene_end), 
                                         abs(temp_gene.gene_start - repeat_elem.repeat_end))

            if temp_distance < distance
                distance = temp_distance
                closest_repeat_start = repeat_elem.repeat_start
                closest_repeat_end = repeat_elem.repeat_end

            end
        end
    end

    te_dist_df[i, :Distance] = distance
    te_dist_df[i, :TE_start] = closest_repeat_start
    te_dist_df[i, :TE_end] = closest_repeat_end
    progress_bits[i] = 1

    if sum(progress_bits) % perc_increment == 0
        println("Progress: ", sum(progress_bits) ÷ perc_increment, "%")

    end
end

# Save the distance df
te_dist_df.Distance = map(dist -> dist == Inf ? missing : Int(dist), te_dist_df.Distance)
te_dist_df.TE_start = map(start -> start == Inf ? missing : Int(start), te_dist_df.TE_start)
te_dist_df.TE_end = map(end_ -> end_ == Inf ? missing : Int(end_), te_dist_df.TE_end)
file_suffix = only_specific_type ? "_$specific_type" : ""
CSV.write("../../dicty_data/te_distance_human$(file_suffix)$(only_h3k9me3 ? "_has_h3k9me3" : "").csv", te_dist_df)
CSV.write("../../dicty_data/te_distance_human_numbered_chroms2$(file_suffix)$(only_h3k9me3 ? "_has_h3k9me3" : "").csv", filter(row -> row.Scaffold in [["$i" for i in 1:22]; "X"; "Y"], te_dist_df))
open("../../dicty_data/cds_ids_num_chroms.txt", "w") do io
    println(io, join(te_dist_df.GeneID, "\n"))
end