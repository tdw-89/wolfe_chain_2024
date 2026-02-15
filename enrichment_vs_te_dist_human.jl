include("prelude.jl")

using Serialization

using .EnrichmentUtils
using .RepeatUtils

# Collect test results for summary table at end
test_results = DataFrame(
    Test=String[],
    PValue=Float64[],
    TestStatistic=Float64[],
    StatisticSymbol=String[],
    Description=String[]
)

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

# Dosage sensitivity data:
significant_overlap_file = "../../dicty_data/lift_over_data_filtered.csv"
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


# TE distance distribution among filtered coding genes vs paralogs without H3K9me3:
paralog_ids = vcat(paralog_data.GeneID, paralog_data.ParalogID)
te_distances = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∉ paralog_ids && id ∉ ids_with_k9me3, te_dist_df.GeneID)]
te_distances_dups = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∈ paralog_ids && id ∉ ids_with_k9me3, te_dist_df.GeneID)]

# What proportion of non-paralog genes overlap a TE vs. paralog genes?
cont_table = [count(d -> d == 0, te_distances) count(d -> d == 0, te_distances_dups);
              count(d -> d > 0, te_distances) count(d -> d > 0, te_distances_dups)]

# Are the two distributions significantly different? (A: Yes - All TEs) <- 0.0010287722758197094
record_test!(test_results,
    "MannWhitneyUTest (TE distance; non-paralog vs paralog; excluding H3K9me3)",
    MannWhitneyUTest(te_distances, te_distances_dups),
    "Mann–Whitney U test: compares TE-distance distributions between non-paralog and paralog genes (excluding genes with H3K9me3); U is the rank-sum-based test statistic.")

# Are the proportions of genes overlapping a TE significantly different between paralogs and non-paralogs? (A: Yes - All TEs) -> 0.0021563409800988667
record_test!(test_results,
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
record_test!(test_results,
    "MannWhitneyUTest (TE distance; with vs without H3K9me3; non-paralogs)",
    MannWhitneyUTest(te_distances_k9me3, te_distances_no_k9me3),
    "Mann–Whitney U test: compares TE-distance distributions between non-paralog genes with H3K9me3 vs without.")
cont_table = [count(d -> d == 0, te_distances_k9me3) count(d -> d == 0, te_distances_no_k9me3); 
              count(d -> d > 0, te_distances_k9me3) count(d -> d > 0, te_distances_no_k9me3)]

# Are the proportions of all genes overlapping a TE significantly different between genes with and without K9me3? (A: Yes - All TEs) -> 0.0021563409800988667
record_test!(test_results,
    "FisherExactTest (TE overlap vs no overlap; with vs without H3K9me3; non-paralogs)",
    FisherExactTest(cont_table[1, 1], cont_table[1, 2], cont_table[2, 1], cont_table[2, 2]),
    "Fisher's exact test on a 2×2 table: tests whether the TE-overlap proportion (Distance==0) differs between non-paralog genes with H3K9me3 vs without.")

# Are the distributions of TE distances significantly different between paralogs with and without K9me3? (A: Yes - All TEs) -> 1.9360047325770461e-59
te_distances_k9me3_dups = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∈ ids_with_k9me3 && id ∈ paralog_ids, te_dist_df.GeneID)]
te_distances_no_k9me3_dups = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> !(id ∈ ids_with_k9me3) && id ∈ paralog_ids, te_dist_df.GeneID)]
fig = plot([box(y=log10.(te_distances_k9me3_dups .+ 1), name="TE Distance w/ K9me3 (Paralogs)", marker_color="blue"), 
      box(y=log10.(te_distances_no_k9me3_dups .+ 1), name="TE Distance w/o k9me3 (Paralogs)", marker_color="red")],
      Layout(yaxis=attr(title="log10(TE Distance)"), plot_bgcolor="rgba(0,0,0,0)"))
record_test!(test_results,
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

record_test!(test_results,
    "FisherExactTest (TE overlap vs no overlap; paralogs with vs without H3K9me3)",
    FisherExactTest(cont_table[1, 1], cont_table[1, 2], cont_table[2, 1], cont_table[2, 2]),
    "Fisher's exact test on a 2×2 table: tests whether the TE-overlap proportion (Distance==0) differs between paralogs with H3K9me3 vs without.")

# Triplosensitivity vs. TE distance vs. K9me3:
significant_overlap_df = CSV.read(significant_overlap_file, DataFrame)
select!(significant_overlap_df, ["gene_id_GRCh38", "gene_name_liftover_GRCh37", "pTriplo_liftover_GRCh37"])
rename!(significant_overlap_df, Dict(:gene_id_GRCh38 => :GeneID, :gene_name_liftover_GRCh37 => :Gene, :pTriplo_liftover_GRCh37 => :pTriplo))

# Individual gene sensitivity vs. TE distance vs. K9me3:
indiv_df = DataFrame(
    :GeneID => vcat(paralog_data.GeneID, paralog_data.ParalogID)
)
indiv_df = innerjoin(indiv_df, significant_overlap_df, on=:GeneID)

insertcols!(indiv_df, :TEDistance => 0.0, :HasK9me3 => false)
indiv_df.TEDistance = [only(te_dist_df.Distance[te_dist_df.GeneID .== gene_id]) for gene_id in indiv_df.GeneID]
indiv_df.HasK9me3 = [gene_id ∈ ids_with_k9me3 for gene_id in indiv_df.GeneID]

display(plot([
    box(y=indiv_df.pTriplo[indiv_df.HasK9me3 .== true], name="With K9me3", marker_color="blue"),
    box(y=indiv_df.pTriplo[indiv_df.HasK9me3 .== false], name="Without K9me3", marker_color="red")
],
Layout(plot_bgcolor="rgba(0,0,0,0)",
yaxis=attr(gridcolor="lightgray",
           gridwidth=1.5))
))

# Are paralogs with H3K9me3 more likely to be triplosensitive? (A: No - 0.6909564111737734)
record_test!(test_results,
    "MannWhitneyUTest (pTriplo; paralogs with vs without H3K9me3)",
    MannWhitneyUTest(indiv_df.pTriplo[indiv_df.HasK9me3 .== true], indiv_df.pTriplo[indiv_df.HasK9me3 .== false]),
    "Mann–Whitney U test: compares triplosensitivity scores (pTriplo) between paralogs with H3K9me3 vs without.")

display(plot([
    box(y=indiv_df.pTriplo[indiv_df.TEDistance .<= 0], name="TE overlap", marker_color="blue"),
    box(y=indiv_df.pTriplo[indiv_df.TEDistance .> 0], name="No TE overlap", marker_color="red")
],
Layout(plot_bgcolor="rgba(0,0,0,0)",
yaxis=attr(gridcolor="lightgray",
           gridwidth=1.5))))

# Are genes with TE overlap more likely to be triplosensitive? (A: Yes - All TEs -> 8.059948448185499e-6)
record_test!(test_results,
    "MannWhitneyUTest (pTriplo; TE overlap vs no TE overlap)",
    MannWhitneyUTest(indiv_df.pTriplo[indiv_df.TEDistance .<= 0], indiv_df.pTriplo[indiv_df.TEDistance .> 0]),
    "Mann–Whitney U test: compares pTriplo between genes overlapping a TE (TEDistance<=0) vs not overlapping (TEDistance>0).")

cont_table = [count(sens -> sens >= 0.5, indiv_df.pTriplo[indiv_df.TEDistance .<= 0]) count(sens -> sens < 0.5, indiv_df.pTriplo[indiv_df.TEDistance .<= 0]); 
              count(sens -> sens >= 0.5, indiv_df.pTriplo[indiv_df.TEDistance .> 0]) count(sens -> sens < 0.5, indiv_df.pTriplo[indiv_df.TEDistance .> 0])]

record_test!(test_results,
    "FisherExactTest (pTriplo≥0.5 vs <0.5; TE overlap vs no overlap)",
    FisherExactTest(cont_table[1, 1], cont_table[1, 2], cont_table[2, 1], cont_table[2, 2]),
    "Fisher's exact test on a 2×2 table: tests whether the fraction of triplosensitive genes (pTriplo≥0.5) differs between TE-overlapping vs non-overlapping genes.")

# --- Summary table of hypothesis tests ---
display(test_results)

out_path = joinpath(@__DIR__, "../../dicty_data", "test_results_enrichment_vs_te_dist_human.csv")
CSV.write(out_path, test_results)