include("prelude.jl")

using .EnrichmentUtils
using .RepeatUtils
using CategoricalArrays

# Collect test results for summary table at end
test_results = DataFrame(
    Test=String[],
    PValue=Float64[],
    TestStatistic=Float64[],
    StatisticSymbol=String[],
    Description=String[]
)

function get_effect_size(test)
    if hasproperty(test, :U)
        U = test.U
        Z = (U - test.mu) / test.sigma |> abs
        N = test.ny + test.nx
        r = Z / √N
        return (r, "r")
    elseif hasproperty(test, :ω)
        return (Float64(getproperty(test, :ω)), "ω")  # odds ratio
    else
        error("Unsupported test type for statistic extraction: $(typeof(test))")
    end
end

function record_test!(results::DataFrame, label::AbstractString, test, description::AbstractString)
    stat_val, stat_sym = get_effect_size(test)
    push!(results, (String(label), Float64(pvalue(test)), stat_val, stat_sym, String(description)))
    return test
end

# Duplicate filtering parameter:
ds_threshold = 3

# Peak files
# chip_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_1_ensembl52/"
chip_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_1_ensembl52/"
atac_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_2_ensembl52/"

# Genome data
gff_source = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
paralog_file = "../../dicty_data/filtered/paralog_filt.tsv"
final_gene_list = "../../dicty_data/filtered/final_gene_list.txt"
te_dist_file = "../../dicty_data/te_distance.csv"
te_density_file = "../../dicty_data/te_density/te_density.csv"
te_id_file = "blacklists/ensembl_te_ids_with_predicted.tsv"

font_family = "Times New Roman"
box_plt_blue(data, name) = box(y=data, 
                        name=name, 
                        marker_color="navy", 
                        line=attr(color="navy"), 
                        fillcolor="#596fff")
box_plt_red(data, name) = box(y=data, 
                        name=name, 
                        marker_color="red", 
                        line=attr(color="red"), 
                        fillcolor="#fc8b8b")
box_layout = Layout(titlefont=attr(family=font_family, size=40),
                    plot_bgcolor="white",
                    yaxis=attr(title="TE Density",
                                titlefont=attr(family=font_family, size=30),
                                tickfont=attr(family=font_family, size=30),
                                gridcolor="lightgray",
                                gridwidth=1.5), 
                    xaxis=attr(titlefont=attr(family=font_family, size=30),
                                tickfont=attr(family=font_family, size=30),
                                gridcolor="lightgray",
                                gridwidth=1.5), 
                    showlegend=false)

# Load the genome data:
ref_genome = loadgenome(gff_source, chrom_lengths_file)

# Load the TE density data:
te_density_df = CSV.read(te_density_file, DataFrame)
rename!(te_density_df, :Gene_ID => :GeneID)
te_density_df.TotalDensity = te_density_df.Upstream .+ te_density_df.Internal .+ te_density_df.Downstream

# Load TE IDs
te_ids = CSV.read(te_id_file, DataFrame)

# Load the paralog data:
paralog_data = CSV.read(paralog_file, DataFrame)

# Load final filtered gene list:
final_gene_list = open(final_gene_list) do file
    readlines(file)
end

# Filter the paralogs:
select!(paralog_data, ["GeneID", "ParalogID", "dS"])
filter!(row -> row["dS"] <= ds_threshold, paralog_data)

# Load peak data
peak_files = [readdir(chip_peak_file_dir, join=true); readdir(atac_peak_file_dir, join=true)]
peak_files = filter(fn -> endswith(fn, ".narrowPeak"), peak_files)
peak_files = filter(fn -> !contains(fn, r"_S[AB]+_"), peak_files)
peak_data = binpeaks(peak_files, chrom_lengths_file)
addtogenes!(ref_genome, peak_data)

# Which genes have H3K9me3?
ids_with_k9me3 = []
for id in te_density_df.GeneID
    gene = get(ref_genome, id)
    
    hasK9me3 = any([
        sum(getsiginrange(gene, GeneRange(TSS(), TES(), -500, 500), ind)) > 0 for ind in [3, 6, 9]
    ])
    if hasK9me3
        push!(ids_with_k9me3, id)
    end
end

# TE density distribution among filtered coding genes vs paralogs:
paralog_ids = vcat(paralog_data.GeneID, paralog_data.ParalogID)
low_ds_paralogs = paralog_data[levelcode.(cut(paralog_data.dS, 10)) .<= 2,:]
paralog_ids_low_ds = vcat(low_ds_paralogs.GeneID, low_ds_paralogs.ParalogID)
te_densities = te_density_df.TotalDensity[map(id -> id ∉ paralog_ids && id ∉ ids_with_k9me3, te_density_df.GeneID)]
te_densities_dups = te_density_df.TotalDensity[map(id -> id ∈ paralog_ids && id ∉ ids_with_k9me3, te_density_df.GeneID)]
fig = plot([box_plt_blue(te_densities, "TE Density"), box_plt_red(te_densities_dups, "TE Density (Paralogs)")],
      box_layout)
display(fig)
savefig(fig, joinpath(@__DIR__, "../../dicty_data", "te_density_nonparalog_vs_paralog.html"))

# Are the TE density distributions significantly different between non-paralogs and paralogs? (A: Yes, p value <1e-06)
record_test!(test_results,
    "MannWhitneyUTest (TE density; non-paralog vs paralog)",
    MannWhitneyUTest(te_densities, te_densities_dups),
    "Mann–Whitney U test: compares TE-density distributions between non-paralog and paralog genes.")

# What is the distribution of TE densities by H3K9me3 status?
te_densities_k9me3 = te_density_df.TotalDensity[map(id -> id ∈ ids_with_k9me3 && id ∉ paralog_ids, te_density_df.GeneID)]
te_densities_no_k9me3 = te_density_df.TotalDensity[map(id -> id ∉ ids_with_k9me3 && id ∉ paralog_ids, te_density_df.GeneID)]

fig = plot([box_plt_blue(te_densities_k9me3, "With H3K9me3"), box_plt_red(te_densities_no_k9me3, "Without H3K9me3")], box_layout)
display(fig)
savefig(fig, joinpath(@__DIR__, "../../dicty_data", "te_density_by_k9me3_status.html"))
record_test!(test_results,
    "MannWhitneyUTest (TE density; non-paralogs with vs without H3K9me3)",
    MannWhitneyUTest(te_densities_k9me3, te_densities_no_k9me3),
    "Mann–Whitney U test on TE density: compares non-paralog genes with H3K9me3 vs without.")

# What is the distribution of distance to the nearest TE


te_densities_k9me3_dups = te_density_df.TotalDensity[map(id -> id ∈ ids_with_k9me3 && id ∈ paralog_ids, te_density_df.GeneID)]
te_densities_no_k9me3_dups = te_density_df.TotalDensity[map(id -> !(id ∈ ids_with_k9me3) && id ∈ paralog_ids, te_density_df.GeneID)]

fig = plot([box_plt_blue(te_densities_k9me3_dups, "With H3K9me3 (Paralogs)"), box_plt_red(te_densities_no_k9me3_dups, "Without H3K9me3 (Paralogs)")], box_layout)
display(fig)
savefig(fig, joinpath(@__DIR__, "../../dicty_data", "te_density_paralogs_by_k9me3_status.html"))
record_test!(test_results,
    "MannWhitneyUTest (TE density; paralogs with vs without H3K9me3)",
    MannWhitneyUTest(te_densities_k9me3_dups, te_densities_no_k9me3_dups),
    "Mann–Whitney U test on TE density: compares paralogs with H3K9me3 vs without.")

# Low dS paralogs
te_densities_low_ds = te_density_df.TotalDensity[map(id -> id ∈ paralog_ids_low_ds, te_density_df.GeneID)]
te_densities_low_ds_k9me3 = te_density_df.TotalDensity[map(id -> id ∈ paralog_ids_low_ds && id ∈ ids_with_k9me3, te_density_df.GeneID)]
te_densities_low_ds_no_k9me3 = te_density_df.TotalDensity[map(id -> id ∈ paralog_ids_low_ds && id ∉ ids_with_k9me3, te_density_df.GeneID)]

fig = plot([box_plt_blue(te_densities_low_ds_k9me3, "With H3K9me3"), box_plt_red(te_densities_low_ds_no_k9me3, "Without H3K9me3")], box_layout)
display(fig)
savefig(fig, joinpath(@__DIR__, "../../dicty_data", "te_density_low_ds_paralogs_by_k9me3_status.html"))
record_test!(test_results,
    "MannWhitneyUTest (TE density; low ds paralogs with vs without H3K9me3)",
    MannWhitneyUTest(te_densities_low_ds_k9me3, te_densities_low_ds_no_k9me3),
    "Mann–Whitney U test on TE density: compares low ds paralogs with H3K9me3 vs without."
)

# --- Summary table of hypothesis tests ---
display(test_results)

out_path = joinpath(@__DIR__, "../../dicty_data", "test_results_enrichment_vs_te_density.csv")
CSV.write(out_path, test_results)