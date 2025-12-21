include("prelude.jl")

using FastaIO
using Serialization

# custom lib:
include("custom_lib/load_gff.jl")
include("custom_lib/enrichment_utils.jl")
include("custom_lib/te_utils.jl")
include("custom_lib/misc_utils.jl")

reload_peak_data = false

# Files
human_gff = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.99.gff3"
chrom_lengths_file = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/chromosome_lengths.txt"
repeat_file = "../../dicty_data/mammals/primates/h_sapiens/GRCh38_UCSC_rmsk.tsv.gz" # UCSC-derived table of repeatmasker-identified repeats
ensembl_cds_file = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.cds.all.fa.gz"
peak_data_dir = "../../dicty_data/mammals/primates/h_sapiens/ENCODE_histone_mods/"
chrom_lengths_file = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/chromosome_lengths.txt"

# Add the peak data
if reload_peak_data
    peak_files = reduce(vcat, [map(fn -> joinpath(root, fn), files) for (root, dir, files) in walkdir(chip_peak_file_dir)])
    peak_files = filter(fn -> endswith(fn, ".bed") || endswith(fn, ".bed.gz"), peak_files)
    peak_files = filter(fn -> contains(fn, "k9me3"), peak_files)
    peak_data = binpeaks(peak_files, chrom_lengths_file)
    serialize("../../dicty_data/julia_serialized/human_h3k9me3_exper.jls", peak_data)
else
    peak_data = deserialize("../../dicty_data/julia_serialized/human_h3k9me3_exper.jls")
end

# Average the peak data
peak_data = average_peak_replicate_groups(peak_data)

# Are there any peaks on non-numbered chromosomes? A: No
main_chroms = [["$i" for i in 1:22]..., "X", "Y"]
peak_data.samples[1].chroms = filter(chr -> chr.name in main_chroms, peak_data.samples[1].chroms)
sample_chrom_names = [chr.name for chr in peak_data.samples[1].chroms]

# Load the repeat df
te_df = CSV.read(repeat_file, DataFrame)

# Filter out low-confidence classes
all_classes = filter(class -> !contains(class, "?"), unique(te_df.repClass))
filter!(row -> row.repClass in all_classes, te_df)

# Select and rename
select!(te_df, [6, 7, 8, 12, 13])
rename!(te_df, [:Chromosome, :Start, :End, :Type, :Family])
te_df.Start = te_df.Start .+ 1 # Convert to 1-based indexing
te_df.Chromosome = map(chr -> String(split(chr, "chr")[2]), te_df.Chromosome)

# Filter to only include the main chromosomes
filter!(row -> row.Chromosome in main_chroms, te_df)

# Calculate the % of each TE family overlapping an h3k9me3 peak
class_family_pairs = unique(zip(te_df.Type, te_df.Family))
coverage_df_family = DataFrame(
                        "Type" => [pair[1] for pair in class_family_pairs],
                        "Family" => [pair[2] for pair in class_family_pairs],
                        "Coverage" => zeros(length(class_family_pairs)),
                        "Overlap" => zeros(length(class_family_pairs)),
                        # "TotalBases" => zeros(Int, length(all_classes)),
                        "TotalN" => zeros(Int, length(class_family_pairs)),
                        "TotalBases" => zeros(Int, length(class_family_pairs))
)

Threads.@threads for i in eachindex(class_family_pairs)
    family = class_family_pairs[i][2]
    family_inds = findall(te_df.Family .== family)
    family_df = te_df[family_inds, :]
    total_n = nrow(family_df)
    coverage_df_family.TotalN[i] = total_n
    overlap_count = 0
    total_signal_sum = 0
    total_bases = sum(family_df.End .- family_df.Start .+ 1)
    coverage_df_family.TotalBases[i] = total_bases

    for (ind, (te_chrom, te_start, te_end, te_type)) in enumerate(zip(family_df.Chromosome, family_df.Start, family_df.End, family_df.Type))
        has_k9me3 = false
        chrom_ind = findfirst(sample_chrom_names .== te_chrom)
        chrom = peak_data.samples[1].chroms[chrom_ind]

        if te_end > length(chrom.signal) # Truncate the TE if it extends past the Y chromosome end
            signal_sum = sum(view(chrom.signal, te_start:length(chrom.signal)))
            has_k9me3 = signal_sum > 0
            overlap_count += has_k9me3
            total_signal_sum += signal_sum
            
        else
            signal_sum = sum(view(chrom.signal, te_start:te_end))
            has_k9me3 = signal_sum > 0
            overlap_count += has_k9me3
            total_signal_sum += signal_sum

        end
    end
    
    overlap_avg = overlap_count / total_n
    coverage_df_family.Overlap[i] = overlap_avg
    coverage_df_family.Coverage[i] = total_signal_sum / total_bases

end

sort!(coverage_df_family, :Coverage, rev=true)
serialize("../../dicty_data/julia_serialized/coverage_df3.jls", coverage_df_family)

# Caclulate overlap proportion for coding genes:

# Load the reference
human_ref = loadgenome(human_gff, chrom_lengths_file)

# Get the CDS IDs
cds_ids = String[]
FastaReader(ensembl_cds_file) do fr
    for (name, seq) in fr
        if occursin("gene_biotype:protein_coding", name) # Excludes T-cell receptor genes, Imunnoglobin genes, and 'processed_pseudogenes'
            push!(cds_ids, match(r"gene:[A-Za-z0-9]+", name).match[6:end])
        
        end
    end
end

# Collect the CDS genes from the reference object
cds_genes = collect(skipmissing(get(human_ref, unique(cds_ids))))
cds_genes = filter(gene -> gene.scaffold.name ∈ main_chroms, cds_genes)

# Calculate the % of all coding sequences covered by h3k9me3 peaks
push!(coverage_df_family, ["CDS", "CDS", 0.0, 0.0, 0, 0])
overlap_count = 0
total_signal_sum = 0
total_bases = sum([gene.gene_end - gene.gene_start + 1 for gene in cds_genes])
coverage_df_family.TotalBases[end] = total_bases
coverage_df_family.TotalN[end] = length(cds_genes)
total_signal_sum = 0

for gene in cds_genes
    gene_start = gene.gene_start
    gene_end = gene.gene_end
    gene_chrom = gene.scaffold.name
    chrom_ind = findfirst(sample_chrom_names .== gene_chrom)
    chrom = peak_data.samples[1].chroms[chrom_ind]
    
    if gene_end > length(chrom.signal) # Truncate the TE if it extends past the Y chromosome end
        signal_sum = sum(view(chrom.signal, gene_start:length(chrom.signal)))
        
    else
        signal_sum = sum(view(chrom.signal, gene_start:gene_end))
        
    end

    overlap_count += signal_sum > 0
    total_signal_sum += signal_sum

end

coverage_df_family.Coverage[end] = total_signal_sum / coverage_df_family.TotalBases[end]
coverage_df_family.Overlap[end] = overlap_count / coverage_df_family.TotalN[end]
# sort!(coverage_df_family, :Coverage, rev=true)
CSV.write("../../dicty_data/h3k9me3_coverage_all.csv", coverage_df_family)

coverage_gdf_family = groupby(coverage_df_family, :Type)
coverage_df_type = combine(coverage_gdf_family, [:Coverage, :Overlap, :TotalN, :TotalBases] => ((c, o, n, b) -> (Coverage=mean(c, Weights(n)), Overlap=mean(o, Weights(b)), TotalN=sum(n), TotalBases=sum(b)) ) => AsTable)
sort!(coverage_df_type, :Coverage, rev=true)
CSV.write("../../dicty_data/h3k9me3_coverage_type.csv", coverage_df_type)