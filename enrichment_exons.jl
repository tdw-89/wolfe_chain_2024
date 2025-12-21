# NOTE: The 'intersect' method from the 'IntervalTrees' package performs redundant comparisons, 
# so there may be a way to speed it up by 2x. For each test of overlep Subject->Query, it also
# tests Query->Subject. It doesn't change the time complexity, but it does double the number of
# comparisons.
include("prelude.jl")

using Serialization
using IntervalTrees
using Graphs
using MetaGraphs

include("custom_lib/genomic_data.jl")

const GENE_TYPE = "singletons"
const RUN_COMPARISON = false

###############
#= FUNCTIONS =#
###############

function filter_to_coding(transcript_df, coding_ids)

    id_strs = [match(r"gene:[A-Z0-9]+", attr).match[6:end] for attr in transcript_df.Attributes]
    contains_coding_id = zeros(Bool, length(id_strs))

    for (i,id) in enumerate(id_strs)

        if id in coding_ids
            contains_coding_id[i] = true
        
        end
    end

    return transcript_df[contains_coding_id, :]
end

function filter_to_supported(transcript_df, support_level)
    if support_level == 1
        support_levels = occursin.(r"transcript_support_level=1", transcript_df.Attributes)

    else
        support_levels = occursin.(Regex("transcript_support_level=[1-$support_level]"), transcript_df.Attributes)

    end

    return transcript_df[support_levels,:]
end

function trim_to_genes(transcript_df, gene_ids)
    
    gene_strs = map(attr_str -> match(r"gene:[A-Z0-9]+", attr_str), transcript_df.Attributes)
    attrs_with_genes = findall(mt -> !isnothing(mt), gene_strs)
    new_df = transcript_df[attrs_with_genes, :]
    gene_id_strs = [attr.match[6:end] for attr in filter(attr -> !isnothing(attr), gene_strs)]
    insertcols!(new_df, :GeneID => gene_id_strs)

    return filter!(row -> row.GeneID in gene_ids, new_df)
end

function merge_intervals(pair::Tuple{IV, IV}) where IV <: IntervalValue
    return IntervalValue(min(first(pair[1]), first(pair[2])), max(last(pair[1]), last(pair[2])), 0)

end

function merge_intervals(interval_1::A, interval_2::A) where A <: AbstractInterval
    return Interval(min(first(interval_1), first(interval_2)), max(last(interval_1), last(interval_2)))

end


"""
    merge_overlaps!(tree::IntervalTree)

Merge overlapping intervals in the tree, and replace them with a single interval that spans the range of the overlapping intervals.
"""
function merge_overlaps!(tree::IntervalTree)
    overlap_itr = intersect(tree, tree)
    
    # Half of the intersections are redundant, so filter to just the pairs where the first interval's value is less than the second's.
    overlaps = filter(pair -> value(pair[1]) < value(pair[2]), collect(overlap_itr))

    # Build a graph of overlapping intervals
    graph = SimpleGraph(length(overlaps))
    meta_graph = MetaGraph(graph)

    for overlap in overlaps
        add_edge!(meta_graph, value(overlap[1]), value(overlap[2]))

        # Will be slightly redundant, but it's probably faster than getting unique values/intervals first and setting them all at once? 
        set_prop!(meta_graph, value(overlap[1]), :interval, Interval(overlap[1]))
        set_prop!(meta_graph, value(overlap[2]), :interval, Interval(overlap[2]))
    
    end

    # Find the connected components
    components = connected_components(meta_graph)

    # Merge the intervals in each component
    for component in components

        if length(component) == 1
            continue

        end

        component_intervals = [get_prop(meta_graph, vertex, :interval) for vertex in component]
        merged_interval = IntervalValue(reduce(merge_intervals, component_intervals), 0)

        for interval in component_intervals
            delete!(tree, interval)
        
        end

        push!(tree, merged_interval)

    end
end

# IntervalTrees.jl overloads:
IntervalTrees.Interval(pair::Tuple{K,K}) where K = Interval(K(pair[1]), K(pair[2]))
IntervalTrees.Interval(interval::IntervalValue) = Interval(first(interval), last(interval))
IntervalTrees.IntervalValue(interval::Interval, value) = IntervalValue(first(interval), last(interval), value)
Base.delete!(tree::IntervalTree, interval::IntervalValue) = Base.delete!(tree, Interval(interval))

# Files
peak_data_serialized_file = "../../dicty_data/julia_serialized/human_h3k9me3_exper.jls"
num_chrom_cds_id_file = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/cds_ids_num_chroms.txt"
transcript_gtf_file = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.99_transcripts.bed"
paralog_file = "../../dicty_data/filtered/human_paralog_info_filt.csv"
singleton_file = "../../dicty_data/filtered/human_singletons_filt.csv"

# Load the numbered chromosome coding gene IDs
coding_ids = readlines(num_chrom_cds_id_file)

# Load the peak data
peak_data = deserialize(peak_data_serialized_file)

# Average the peak data
peak_data = average_peak_replicate_groups(peak_data)
main_chroms = [["$i" for i in 1:22]..., "X", "Y"]
peak_data.samples[1].chroms = filter(chr -> chr.name in main_chroms, peak_data.samples[1].chroms)

# Load the transcript exon data
transcript_df = CSV.read(transcript_gtf_file, DataFrame, header=false, delim="\t")
rename!(transcript_df, [:Chrom, 
                        :TranscriptStart, 
                        :TranscriptEnd, 
                        :TranscriptID, 
                        :Score, 
                        :Strand, 
                        :ThickStart,
                        :ThickEnd,
                        :RGB,
                        :nExon,
                        :Widths,
                        :Starts,
                        :Attributes])

# Filter to coding genes on numbered chromosomes consistent with those used for TE distances
transcript_df = filter_to_coding(transcript_df, coding_ids)

# Filter to transcript support level <= 2
transcript_df = filter_to_supported(transcript_df, 2)

# If ONLY_DUPS is true, filter to only the duplicates
if GENE_TYPE == "duplicates"
    paralog_df = CSV.read(paralog_file, DataFrame)
    paralog_ids = unique(vcat(paralog_df.GeneID, paralog_df.ParalogID))
    transcript_df = trim_to_genes(transcript_df, paralog_ids)

elseif GENE_TYPE == "singletons"
    singleton_df = CSV.read(singleton_file, DataFrame)
    singleton_ids = singleton_df.GeneID
    transcript_df = trim_to_genes(transcript_df, singleton_ids)

end

# Convert the transcript list to a BED3-style table of exons

    # Count the total number of exons, and pre-allocate three vectors
total_exon_n = sum(transcript_df.nExon)
chroms = mapreduce((chrom, exon_n) -> repeat([String(chrom)], exon_n), 
                    vcat, 
                transcript_df.Chrom, 
                transcript_df.nExon)
starts = mapreduce((start, chrom_start) -> 
    parse.(Int, split(start, ",")[1:end-1]) .+ chrom_start .+ 1, 
    vcat, 
    transcript_df.Starts, 
    transcript_df.TranscriptStart)
ends = mapreduce(width -> parse.(Int, split(width, ",")[1:end-1]), 
                vcat, 
                transcript_df.Widths) .+ (starts .- 1)

    # Create the exon table
exon_df = DataFrame(Chromosome = chroms, Start = starts, End = ends)
unique!(exon_df)

    # Sort by chromosome and then start
!issorted(exon_df, [:Chromosome, :Start]) && sort!(exon_df, [:Chromosome, :Start])

tree_dict = Dict(
    begin
    df = exon_df[exon_df.Chromosome .== chr, [:Start, :End]]
    intervals = [IntervalValue{Int, Int}(interval[1], interval[2], i) for (i, interval) in enumerate(zip(df.Start, df.End))]
    chr => IntervalTree{Int, IntervalValue{Int, Int}}(sort(intervals))
    end
    for chr in main_chroms
)

# Merge overlapping exons
for chr in keys(tree_dict)
    merge_overlaps!(tree_dict[chr])

end

# Calculate the % of total exon bases overlapping an H3K9me3 peak in at least one sample:
exon_df_merged = DataFrame()
for (chr,tree) in tree_dict
    df = DataFrame(tree)

    if nrow(df) > 0 
        df.Chromosome = repeat([chr], nrow(df))
        append!(exon_df_merged, df)
    
    end
end

select!(exon_df_merged, [4, 1, 2])
rename!(exon_df_merged, [:Chromosome, :Start, :End])
exon_df_merged.Length = exon_df_merged.End .- exon_df_merged.Start .+ 1
exon_df_merged.OverlapPerc = zeros(Float64, nrow(exon_df_merged))

signal_dict = Dict(chrom.name => chrom.signal for chrom in peak_data.samples[1].chroms)

for i in 1:nrow(exon_df_merged)
    chrom = exon_df_merged.Chromosome[i]
    start = exon_df_merged.Start[i]
    stop = exon_df_merged.End[i]
    chrom_signal = signal_dict[chrom]
    signal_sum = sum(chrom_signal[start:stop])
    exon_df_merged.OverlapPerc[i] = signal_sum / exon_df_merged.Length[i]
    
end

CSV.write("../../dicty_data/human_h3k9me3_exon_enrichment$(GENE_TYPE in ["singletons", "duplicates"] ? "_$GENE_TYPE" : "").csv", exon_df_merged)

# Calculate the weighted mean
total_bases = sum(exon_df_merged.Length)
total_signal_sum = sum(exon_df_merged.OverlapPerc .* exon_df_merged.Length)
weighted_mean = (total_signal_sum / total_bases) * 100

# Comparison of singletons to duplicates:
if RUN_COMPARISON
    using HypothesisTests
    singleton_df = CSV.read("../../dicty_data/human_h3k9me3_exon_enrichment_singletons.csv", DataFrame)
    duplicate_df = CSV.read("../../dicty_data/human_h3k9me3_exon_enrichment_duplicates.csv", DataFrame)
    weighted_vals_singleton = singleton_df.OverlapPerc .* singleton_df.Length
    weighted_vals_duplicate = duplicate_df.OverlapPerc .* duplicate_df.Length
    p_val = pvalue(MannWhitneyUTest(weighted_vals_singleton, weighted_vals_duplicate))

end