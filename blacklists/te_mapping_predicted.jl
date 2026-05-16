#=
STEP 2: Get a list of TE IDs from the DB reference and those overlapping a predicted TE
=#
include("../prelude.jl")

const TE_OVERLAP_CUTOFF = 80

function merge_overlapping_intervals(ranges::Vector{Tuple{Any, Any}})
    if isempty(ranges)
        return Tuple{Int, Int}[]
    end

    # Sort ranges by start position
    sorted_ranges = sort(ranges, by = x -> x[1])
    merged = Tuple{Int, Int}[]
    current_start, current_end = sorted_ranges[1]

    for (range_start, range_end) in sorted_ranges[2:end]
        if range_start <= current_end + 1  # Overlapping or contiguous intervals
            current_end = max(current_end, range_end)
        else
            push!(merged, (current_start, current_end))
            current_start, current_end = range_start, range_end
        end
    end
    push!(merged, (current_start, current_end))  # Add the last interval

    return merged
end


function merge_overlapping_intervals(ranges::Vector{Tuple{Int, Int}})
    merge_overlapping_intervals(Tuple{Any, Any}[(r[1], r[2]) for r in ranges])
end


function overlap_lens(gene::Gene)
    ba(s::Int,e::Int) = zeros(Bool, e - s + 1) |> BitArray
    gene_bit_arr = ba(gene.gene_start, gene.gene_end)
    exon_ranges = [(exon.exon_start, exon.exon_end) for exon in gene.exons]
    exons_ranges_merged = merge_overlapping_intervals(exon_ranges)
    exon_bit_arrs = [ba(exon_range[1], exon_range[2]) for exon_range in exons_ranges_merged]
    for repeat_elem in gene.scaffold.repeats
        if GenomeTypes.hasoverlap(
            repeat_elem.repeat_start,
            gene.gene_start,
            repeat_elem.repeat_end,
            gene.gene_end
        )
            overlap_start = max(repeat_elem.repeat_start, gene.gene_start)
            overlap_end = min(repeat_elem.repeat_end, gene.gene_end)
            for pos in overlap_start:overlap_end
                gene_bit_arr[pos - gene.gene_start + 1] = true
            end

            for (i, exon_range) in enumerate(exons_ranges_merged)
                if GenomeTypes.hasoverlap(
                    repeat_elem.repeat_start,
                    exon_range[1],
                    repeat_elem.repeat_end,
                    exon_range[2]
                )
                    exon_overlap_start = max(repeat_elem.repeat_start, exon_range[1])
                    exon_overlap_end = min(repeat_elem.repeat_end, exon_range[2])
                    for pos in exon_overlap_start:exon_overlap_end
                        exon_bit_arrs[i][pos - exon_range[1] + 1] = true
                    end
                end
            end
        end
    end

    gene_overlap_perc = count(gene_bit_arr) / length(gene_bit_arr) * 100.0
    exon_overlap_percs = sum([count(exon_bit_arr) for exon_bit_arr in exon_bit_arrs]) / sum([length(exon_bit_arr) for exon_bit_arr in exon_bit_arrs]) * 100.0
    return gene_overlap_perc, exon_overlap_percs
end

using CSV
using DataFrames
using SparseArrays
using .RepeatUtils

cd(@__DIR__)

# Genome data:
gff_source = "../../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"

# Predicted TE files (Ensembl 52):
tes_ensembl = "../../../dicty_data/rm_genomes/ensembl_52/full/results_v3/onecodetofindthemall/Dictyostelium_discoideum.dicty_2.7.dna.toplevel_aggregated_compat.csv"

# mapped TEs file:
mapped_tes_file = "ensembl_te_ids_full.tsv"

# Load the genome data:
ref_genome = loadgenome(gff_source, chrom_lengths_file, feature_type="all")

# Load the TE data:
te_df_ensembl = CSV.read(tes_ensembl, DataFrame)
mapped_tes = CSV.read(mapped_tes_file, DataFrame)

# Add the repeats to the genome:
convert_to_repeats!(
    ref_genome,
    te_df_ensembl
    ; allow_missing_scaffolds=false)

# Get genes that overlap with predicted TEs:
te_genes = [findoverlappinggenes(repeat_elem) for repeat_elem in ref_genome.repeats]
overlapping_ids = [te_genes[i][1] for i in eachindex(te_genes) if !ismissing(te_genes[i])]
overlapping_ids = reduce(vcat, overlapping_ids)
overlapping_ids = unique(overlapping_ids)

overlap_perc_df = DataFrame(
    :GeneID => overlapping_ids,
    :TotalOverlapPerc => zeros(Float64, length(overlapping_ids)),
    :ExonOverlapPerc => zeros(Float64, length(overlapping_ids)),
    :IntronOverlapPerc => zeros(Float64, length(overlapping_ids))
)
for (i, gene_id) in enumerate(overlapping_ids)
    gene = ref_genome.genes[2][findfirst(g -> g.id == gene_id, ref_genome.genes[2])]
    total_perc, exon_perc = overlap_lens(gene)
    overlap_perc_df[i, :TotalOverlapPerc] = total_perc
    overlap_perc_df[i, :ExonOverlapPerc] = exon_perc
    overlap_perc_df[i, :IntronOverlapPerc] = max(total_perc - (total_perc * exon_perc / 100), 0.0)
end

sort!(overlap_perc_df, :ExonOverlapPerc)
CSV.write("../../../dicty_data/te_overlap_perc_dicty.csv", overlap_perc_df)
overlapping_ids_filt = filter(row -> row.ExonOverlapPerc >= TE_OVERLAP_CUTOFF, overlap_perc_df).GeneID

# Merge the overlapping IDs with the TE IDs lifted from the dictybase genome:
all_te_IDs = union(mapped_tes.GeneID, overlapping_ids_filt)
CSV.write("ensembl_te_ids_with_predicted.tsv", DataFrame(GeneID=all_te_IDs))
labeled_df = DataFrame(GeneID=mapped_tes.GeneID, Source="dictyBase")
labeled_df = vcat(labeled_df, DataFrame(GeneID=overlapping_ids_filt, Source="Predicted"))
CSV.write("ensembl_te_ids_with_source.tsv", labeled_df)

# Total genome coverage
#
#   Combined: (4205255 / 34134454) = 12.32
#   Class I: (3746292 / 34134454) = 10.97
#   Class II: (458963 / 34134454) = 1.34
for te_type in ["All", "I", "II"]
    te_df_ensembl_copy = deepcopy(te_df_ensembl)
    if te_type == "I"
        filter!(row -> !contains(row.Type, "DNA"), te_df_ensembl_copy)
    elseif te_type == "II"
        filter!(row -> contains(row.Type, "DNA"), te_df_ensembl_copy)
    end
    genome_size = values(ref_genome.scaffolds) |>
        collect |> 
        sv -> map(s -> s.scaffold_end - s.scaffold_start + 1, sv) |> 
        sum

    sum([te.End - te.Start + 1 for te in eachrow(te_df_ensembl_copy)]) / genome_size * 100.0 |> println
end