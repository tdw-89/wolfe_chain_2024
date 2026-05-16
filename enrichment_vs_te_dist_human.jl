include("prelude.jl")

using GaussianMixtures
using MultipleTesting
using Random
using Serialization

using .EnrichmentUtils
using .RepeatUtils

function perm_test(
    stat_func::F,
    group1::Vector{M},
    group2::Vector{N}; n_permutations=10_000
    ) where { F <: Function, M <: Number, N <: Number }
    group2_copy = copy(group2) # To avoid modifying the original group2 during shuffling
    observed_stat = stat_func(group1, group2)
    count_extreme = 0
    for _ in 1:n_permutations
        shuffle!(group2_copy)
        perm_stat = stat_func(group1, group2_copy)
        if abs(perm_stat) >= abs(observed_stat)
            count_extreme += 1
        end
    end
    p_value = (count_extreme + 1) / (n_permutations + 1) # Add 1 to numerator and denominator for continuity correction
    return observed_stat, p_value
end

function get_statistic(test)
    if hasproperty(test, :U)
        return (Float64(getproperty(test, :U)), "U")
    elseif hasproperty(test, :ω)
        return (Float64(getproperty(test, :ω)), "ω")  # odds ratio parameter
    else
        error("Unsupported test type for statistic extraction: $(typeof(test))")
    end
end

function record_test!(results::DataFrame, label::AbstractString, test, description::AbstractString)
    stat_val, stat_sym = get_statistic(test)
    push!(results, (String(label), Float64(pvalue(test)), stat_val, stat_sym, String(description)))
    return test
end

reload_peak_data = true
te_type = "TE"

# Peak files
# chip_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_1_ensembl52/"
chip_peak_file_dir = "../../dicty_data/mammals/primates/h_sapiens/ENCODE_histone_mods/"

# Genome data
gff_source = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.99.gff3"
chrom_lengths_file = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/chromosome_lengths.txt"
paralog_file = "../../dicty_data/filtered/human_paralog_info_filt.csv"
te_dist_file = "../../dicty_data/te_distance_human_$te_type.csv"
te_density_file = "../../dicty_data/te_density/te_density_human.csv"

# Dosage sensitivity data:
te_predictions_file = "../../dicty_data/lift_over_data_filtered.csv"
mismatch_file = "../../dicty_data/lift_over_data_mismatches.csv"

# Load the genome data:
ref_genome = loadgenome(gff_source, chrom_lengths_file)

# Load the TE distance data:
te_dist_df = CSV.read(te_dist_file, DataFrame)

# Load the paralog data:
paralog_data = CSV.read(paralog_file, DataFrame)

# Filter the paralogs:
select!(paralog_data, ["GeneID", "ParalogID", "dS"])
filter!(row -> row.dS <= 3, paralog_data)

fig = plot(histogram(x=paralog_data.dS), Layout(title="Distribution of dS for filtered paralogs", xaxis=attr(title="dS")))
savefig(fig, "../../dicty_data/paralog_dS_histogram.html")

# Load peak data
if reload_peak_data
    peak_files = reduce(vcat, [map(fn -> joinpath(root, fn), files) for (root, dir, files) in walkdir(chip_peak_file_dir)])
    peak_files = filter(fn -> endswith(fn, ".bed") || endswith(fn, ".bed.gz"), peak_files)
    peak_files = filter(fn -> contains(fn, "k9me3"), peak_files)
    peak_data = binpeaks(peak_files, chrom_lengths_file)
    serialize("../../dicty_data/julia_serialized/human_h3k9me3_exper.jls", peak_data)
else
    peak_data = deserialize("../../dicty_data/julia_serialized/human_h3k9me3_exper.jls")
end

addtogenes!(ref_genome, peak_data)

# Which genes have K9me3?
ids_with_k9me3 = []

for id in te_dist_df.GeneID

    gene = get(ref_genome, id)
    
    has_k9 = any([
        any(sum(getsiginrange(gene, GeneRange(TSS(), TES(), -500, 500), ind)) > 0) for ind in eachindex(gene.samples)
    ])
    if has_k9
        push!(ids_with_k9me3, id)
    end
end

# Collect test results per group for within-group BH adjustment
group_a_results = DataFrame(
    Test=String[],
    PValue=Float64[],
    TestStatistic=Float64[],
    StatisticSymbol=String[],
    Description=String[]
)
group_b_results = DataFrame(
    Test=String[],
    PValue=Float64[],
    TestStatistic=Float64[],
    StatisticSymbol=String[],
    Description=String[]
)
# BEGIN GROUP A

# TE distance distribution among filtered coding genes vs paralogs without H3K9me3:
paralog_ids = vcat(paralog_data.GeneID, paralog_data.ParalogID)
te_distances = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∉ paralog_ids && id ∉ ids_with_k9me3, te_dist_df.GeneID)]
te_distances_dups = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∈ paralog_ids && id ∉ ids_with_k9me3, te_dist_df.GeneID)]

# What proportion of non-paralog genes overlap a TE vs. paralog genes?
cont_table = [count(d -> d == 0, te_distances) count(d -> d == 0, te_distances_dups);
              count(d -> d > 0, te_distances) count(d -> d > 0, te_distances_dups)]

# Are the two distributions significantly different? (A: Yes - All TEs) <- 0.0010287722758197094
record_test!(group_a_results,
    "MannWhitneyUTest (TE distance; non-paralog vs paralog; excluding H3K9me3)",
    MannWhitneyUTest(te_distances, te_distances_dups),
    "Mann–Whitney U test: compares TE-distance distributions between non-paralog and paralog genes (excluding genes with H3K9me3); U is the rank-sum-based test statistic.")

# Are the proportions of genes overlapping a TE significantly different between paralogs and non-paralogs? (A: Yes - All TEs) -> 0.0021563409800988667
record_test!(group_a_results,
    "FisherExactTest (TE overlap vs no overlap; non-paralog vs paralog; excluding H3K9me3)",
    FisherExactTest(cont_table[1, 1], cont_table[1, 2], cont_table[2, 1], cont_table[2, 2]),
    "Fisher's exact test on a 2×2 table: tests whether the TE-overlap proportion (Distance==0) differs between non-paralog and paralog genes (excluding genes with H3K9me3).")

# What is the distribution of TE distances among coding genes w/ vs. w/o H3K9me3?
# What about among just the filtered paralogs?
te_distances_k9me3 = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∈ ids_with_k9me3 && id ∉ paralog_ids, te_dist_df.GeneID)]
te_distances_no_k9me3 = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∉ ids_with_k9me3 && id ∉ paralog_ids, te_dist_df.GeneID)]

# Is the distribution of TE distances significantly different between genes with and without K9me3? (A: Yes - All TEs) -> 2.6419280524529124e-78
plot([box(y=te_distances_k9me3, name="TE Distance w/ K9me3", marker_color="blue"), 
      box(y=te_distances, name="TE Distance w/o k9me3", marker_color="red")])
record_test!(group_a_results,
    "MannWhitneyUTest (TE distance; with vs without H3K9me3; non-paralogs)",
    MannWhitneyUTest(te_distances_k9me3, te_distances_no_k9me3),
    "Mann–Whitney U test: compares TE-distance distributions between non-paralog genes with H3K9me3 vs without.")
cont_table = [count(d -> d == 0, te_distances_k9me3) count(d -> d == 0, te_distances_no_k9me3); 
              count(d -> d > 0, te_distances_k9me3) count(d -> d > 0, te_distances_no_k9me3)]

# Are the proportions of all genes overlapping a TE significantly different between genes with and without K9me3? (A: Yes - All TEs) -> 0.0021563409800988667
record_test!(group_a_results,
    "FisherExactTest (TE overlap vs no overlap; with vs without H3K9me3; non-paralogs)",
    FisherExactTest(cont_table[1, 1], cont_table[1, 2], cont_table[2, 1], cont_table[2, 2]),
    "Fisher's exact test on a 2×2 table: tests whether the TE-overlap proportion (Distance==0) differs between non-paralog genes with H3K9me3 vs without.")

# Are the distributions of TE distances significantly different between paralogs with and without K9me3? (A: Yes - All TEs) -> 1.9360047325770461e-59
te_distances_k9me3_dups = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∈ ids_with_k9me3 && id ∈ paralog_ids, te_dist_df.GeneID)]
te_distances_no_k9me3_dups = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> !(id ∈ ids_with_k9me3) && id ∈ paralog_ids, te_dist_df.GeneID)]
fig = plot([box(y=log10.(te_distances_k9me3_dups .+ 1), name="TE Distance w/ K9me3 (Paralogs)", marker_color="blue"), 
      box(y=log10.(te_distances_no_k9me3_dups .+ 1), name="TE Distance w/o k9me3 (Paralogs)", marker_color="red")],
      Layout(yaxis=attr(title="log10(TE Distance)"), plot_bgcolor="rgba(0,0,0,0)"))
record_test!(group_a_results,
    "MannWhitneyUTest (TE distance; paralogs with vs without H3K9me3)",
    MannWhitneyUTest(te_distances_k9me3_dups, te_distances_no_k9me3_dups),
    "Mann–Whitney U test: compares TE-distance distributions between paralogs with H3K9me3 vs without.")
savefig(fig, joinpath(@__DIR__, "../../dicty_data", "te_distance_by_k9me3_status_paralogs_human.html"))

# Are the ratios of TE overlap significantly different between paralogs with and without K9me3? (A: Yes - All TEs) -> 25.493113707490779e-66
plot(bar(x=["Paralogs With H3K9me3", "Paralogs Without H3K9me3"], 
         y=[count(d -> d == 0, te_distances_k9me3_dups)/length(te_distances_k9me3_dups), 
            count(d -> d == 0, te_distances_no_k9me3_dups)/length(te_distances_no_k9me3_dups)]),
            Layout(yaxis=attr(title="% Overlapping a TE", range=[0, 1])))
cont_table = [count(d -> d == 0, te_distances_k9me3_dups) count(d -> d == 0, te_distances_no_k9me3_dups); 
              count(d -> d > 0, te_distances_k9me3_dups) count(d -> d > 0, te_distances_no_k9me3_dups)]

record_test!(group_a_results,
    "FisherExactTest (TE overlap vs no overlap; paralogs with vs without H3K9me3)",
    FisherExactTest(cont_table[1, 1], cont_table[1, 2], cont_table[2, 1], cont_table[2, 2]),
    "Fisher's exact test on a 2\u00d72 table: tests whether the TE-overlap proportion (Distance==0) differs between paralogs with H3K9me3 vs without.")

# Benjamini-Hochberg correction across all Group A p-values
group_a_results.PValue = adjust(group_a_results.PValue, BenjaminiHochberg())

# END GROUP A
# BEGIN GROUP B

# Triplosensitivity vs. TE distance vs. K9me3:
mapped_ts_predictions_df = CSV.read(te_predictions_file, DataFrame)
select!(mapped_ts_predictions_df, ["gene_id_GRCh38", "gene_name_liftover_GRCh37", "pTriplo_liftover_GRCh37"])
rename!(mapped_ts_predictions_df, Dict(:gene_id_GRCh38 => :GeneID, :gene_name_liftover_GRCh37 => :Gene, :pTriplo_liftover_GRCh37 => :pTriplo))

# Individual gene sensitivity vs. TE distance vs. K9me3:
indiv_df = DataFrame(
    :GeneID => vcat(paralog_data.GeneID, paralog_data.ParalogID)
)
indiv_df = innerjoin(indiv_df, mapped_ts_predictions_df, on=:GeneID)

insertcols!(indiv_df, :TEDistance => 0.0, :HasK9me3 => false)
indiv_df.TEDistance = [only(te_dist_df.Distance[te_dist_df.GeneID .== gene_id]) for gene_id in indiv_df.GeneID]
indiv_df.HasK9me3 = [gene_id ∈ ids_with_k9me3 for gene_id in indiv_df.GeneID]

# Join TE density data (upstream, internal, downstream scores per gene)
te_density_df = CSV.read(te_density_file, DataFrame)
rename!(te_density_df, :Gene_ID => :GeneID)
indiv_df = innerjoin(indiv_df, te_density_df, on=:GeneID)

display(plot([
    box(y=indiv_df.pTriplo[indiv_df.HasK9me3 .== true], name="With K9me3", marker_color="blue"),
    box(y=indiv_df.pTriplo[indiv_df.HasK9me3 .== false], name="Without K9me3", marker_color="red")
],
Layout(plot_bgcolor="rgba(0,0,0,0)",
yaxis=attr(gridcolor="lightgray",
           gridwidth=1.5))
))

# Are paralogs with H3K9me3 more likely to be triplosensitive? (A: No - 0.6909564111737734)
record_test!(group_b_results,
    "MannWhitneyUTest (pTriplo; paralogs with vs without H3K9me3)",
    MannWhitneyUTest(indiv_df.pTriplo[indiv_df.HasK9me3 .== true], indiv_df.pTriplo[indiv_df.HasK9me3 .== false]),
    "Mann–Whitney U test: compares triplosensitivity scores (pTriplo) between paralogs with H3K9me3 vs without.")

display(plot([
    box(y=indiv_df.pTriplo[indiv_df.Internal .> 0], name="Internal TE density > 0", marker_color="blue"),
    box(y=indiv_df.pTriplo[indiv_df.Internal .<= 0], name="Internal TE density = 0", marker_color="red")
],
Layout(plot_bgcolor="rgba(0,0,0,0)",
yaxis=attr(gridcolor="lightgray",
           gridwidth=1.5))))

# Are genes with internal TE density more likely to be triplosensitive?
record_test!(group_b_results,
    "MannWhitneyUTest (pTriplo; Internal TE density > 0 vs = 0)",
    MannWhitneyUTest(indiv_df.pTriplo[indiv_df.Internal .> 0], indiv_df.pTriplo[indiv_df.Internal .<= 0]),
    "Mann–Whitney U test: compares pTriplo between paralogs with non-zero internal TE density vs zero internal TE density.")

# Permutation tests: Spearman correlation between pTriplo and each TE density component
ptriplo    = collect(Float64, indiv_df.pTriplo)
upstream   = collect(Float64, indiv_df.Upstream)
downstream = collect(Float64, indiv_df.Downstream)

perm_results = Tuple{String, Tuple{Float64, Float64}, String}[]

# Internal: Spearman correlation among genes with non-zero internal TE density
internal   = collect(Float64, indiv_df.Internal)
push!(perm_results, (
    "Internal (density > 0)",
    perm_test(corspearman, ptriplo, internal),
    "Permutation test (Spearman ρ): correlation between pTriplo and internal TE density"
))

# Upstream and Downstream: classify into two density clusters first, then test each cluster
for (region_name, density) in [("Upstream", upstream), ("Downstream", downstream)]
    push!(perm_results, (
        "$region_name",
        perm_test(corspearman, ptriplo, density),
        "Permutation test (Spearman ρ): correlation between pTriplo and $region_name TE density"
    ))
end

# Push permutation results into group_b_results (raw p-values; adjusted below with all Group B tests)
for i in eachindex(perm_results)
    label, (statistic, raw_p), description = perm_results[i]
    push!(group_b_results, (label, raw_p, statistic, "\u03c1", description))
end

# Benjamini-Hochberg correction across all Group B p-values
group_b_results.PValue = adjust(group_b_results.PValue, BenjaminiHochberg())

# END GROUP B

# Combine all test results
test_results = vcat(group_a_results, group_b_results)

# --- Summary table of hypothesis tests ---
for row in eachrow(sort(test_results, :PValue))
    println(row.Description)
    println("---- 𝑝 value: " * string(row.PValue))
    println("---- " * row.StatisticSymbol * ": " * string(row.TestStatistic))
    println("")
end

out_path = joinpath(@__DIR__, "../../dicty_data", "test_results_enrichment_vs_te_dist_human_$te_type.csv")
CSV.write(out_path, test_results)