include("prelude.jl")

using CategoricalArrays
using NaturalSort
using RollingFunctions
using JSON
using Printf

using .EnrichmentUtils
using .MiscUtils

function adjust_expr(gene::Gene, ratio_df::DataFrame, id_col::Int64)
    ratio = only(ratio_df.InputRatio[ratio_df[:,id_col] .== gene.id])
    return log2((only(only(gene.rnas).expression) / ratio) + 0.5)

end

function euclid_dist(u, v)
    return sqrt(sum((u .- v).^2))

end

function ifany(v)
    if !isempty(v.rnas)
        return v.rnas[1].expression[1]
    else
        return Inf
    end
end

NEWDS = use_yn00()  # --yn00: YN00 dS instead of Ensembl Compara dS

# Plotting parameters
font_family = "Times New Roman"
y_axis_fmt_expr = attr(title = "log2(TPM + 0.5)", 
                titlefont=attr(family=font_family, size=30),
                tickfont=attr(family=font_family, size=30),
                gridcolor="lightgray",
                gridwidth=1.5,
                range=[0,12])
y_axis_fmt_diff = attr(title = "Euclidean Distance", 
                titlefont=attr(family=font_family, size=30),
                tickfont=attr(family=font_family, size=30),
                gridcolor="lightgray",
                gridwidth=1.5)
x_axis_fmt = attr(title = "𝑑𝑆 (Decile)", 
                titlefont=attr(family=font_family, size=30),
                tickfont=attr(family=font_family, size=30),
                gridcolor="lightgray",
                gridwidth=1.5)
box_plt(data, name) = box(y=data, 
                            name=name, 
                            marker_color="blue", 
                            line=attr(color="blue"), 
                            fillcolor="#596fff")
box_layout_expr = Layout(titlefont=attr(family=font_family, size=40),
                    plot_bgcolor="white",
                    yaxis=y_axis_fmt_expr, 
                    xaxis=x_axis_fmt, 
                    showlegend=false)
box_layout_diff = Layout(titlefont=attr(family=font_family, size=40),
                    plot_bgcolor="white",
                    yaxis=y_axis_fmt_diff, 
                    xaxis=x_axis_fmt, 
                    showlegend=false)

# Genome data
gff_data = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"

# Expression data
expr_data_file = "../../dicty_data/filtered/expr_data_filt_kallisto_ensembl52_single.tsv"
expr_data_adj = "../../dicty_data/filtered/expr_data_filt_adj.tsv"

# Peak files
chip_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_1_ensembl52/"
atac_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_2_ensembl52/"

# Paralog data
paralog_file = "../../dicty_data/filtered/paralog_filt.tsv"

# Singleton list
singleton_list_file = "../../dicty_data/filtered/singleton_filt.tsv"
# final_gene_list = "../../dicty_data/filtered/final_gene_list.txt"

# Load expression data
expr_data = CSV.read(expr_data_file, DataFrame)

# Rename the expression columns to a letter indicating with life-cycle stage the sample
# was taken from.
rename!(expr_data, ["GeneID", "V", "S", "M", "F"])
select!(expr_data, ["GeneID", "V", "M", "F"])  # stages with ChIP-seq data
expr_data[!,2:end] = log2.(expr_data[!,2:end] .+ 0.5)

# Load peak data
peak_files = [filter(fn -> !contains(fn, r"_S[AB]_"), readdir(chip_peak_file_dir, join=true)); readdir(atac_peak_file_dir, join=true)]
filter!(fn -> endswith(fn, ".narrowPeak"), peak_files)
peak_data = binpeaks(peak_files, chrom_lengths_file)

# Load paralog and singleton data
paralog_data = CSV.read(paralog_file, DataFrame)
if NEWDS
    # The same pairs with YN00 dS (written by ds_dist_comparison.qmd)
    paralog_data = CSV.read("../../dicty_data/old_pairs_with_new_ds_dicty.csv", DataFrame)
    paralog_data.dS = paralog_data.New_dS
    filter!(row -> !ismissing(row["dS"]), paralog_data)
    paralog_data.dS = Float64.(paralog_data.dS)
end
singleton_list = CSV.read(singleton_list_file, DataFrame)

# Filter the singleton data to just those that have expression data (paralogs are already filtered)
filter!(row -> row.GeneID in expr_data.GeneID, singleton_list)

# Filter paralogs according to id threshold
select!(paralog_data, ["GeneID", "ParalogID", "dS"])
filter!(row -> row.dS <= 3, paralog_data)

# Filter to pairs which have expression data for both
filter!(row -> row.GeneID in expr_data.GeneID && row.ParalogID in expr_data.GeneID, paralog_data) 

# Load reference genome object
ref_genome = loadgenome(gff_data, chrom_lengths_file)
ref_genome_adj = loadgenome(gff_data, chrom_lengths_file)

# Add data to reference
addexpression!(ref_genome, expr_data)
addexpression!(ref_genome_adj, filter(row -> !ismissing(row.Avg), CSV.read(expr_data_adj, DataFrame)))
addtogenes!(ref_genome, peak_data)
addtogenes!(ref_genome_adj, peak_data)

# Add pair-average expression to paralog data
insertcols!(paralog_data, :AvgExpr => 
    [mean([mean(get(ref_genome, pair[1]).rnas[1].expression), mean(get(ref_genome, pair[2]).rnas[1].expression)]) for pair in zip(paralog_data.GeneID, paralog_data.ParalogID)]
)

insertcols!(paralog_data, :AvgExprAdj => 
    [mean([mean(ifany(get(ref_genome_adj, pair[1]))), 
           mean(ifany(get(ref_genome_adj, pair[2])))]) for pair in zip(paralog_data.GeneID, 
                                                                       paralog_data.ParalogID)]
)

# Add expression difference to paralog data
insertcols!(paralog_data, :Diff => 
    [euclid_dist(get(ref_genome, pair[1]).rnas[1].expression, get(ref_genome, pair[2]).rnas[1].expression) for pair in zip(paralog_data.GeneID, paralog_data.ParalogID)]
)

sort!(paralog_data, :dS)
cor_test = perm_cor_2side(paralog_data.AvgExpr, Vector(paralog_data.dS))

# Get the dS values and quantiles for the filtered pairs
quantile_labels = cut(paralog_data[!, 3], 10)
# insertcols!(paralog_data, :MaxPerc => mean(hcat(paralog_data[:,3], paralog_data[:,4]), dims=2)[:,1])
# quantile_labels = cut(paralog_data.MaxPerc, 10)
quantile_vals = levelcode.(quantile_labels)

# Get the significance of differences between median pair-average expression in each decile and median singleton expression
singletons =  get(ref_genome, Vector(singleton_list.GeneID))
# Stage-averaged, like the paralog AvgExpr
singleton_expr_vals = [mean(gene.rnas[1].expression) for gene in singletons if !ismissing(gene) && !isempty(gene.rnas)]
expression_deciles = [paralog_data.AvgExpr[quantile_vals .== q] for q in sort(unique(quantile_vals))]

open("../../dicty_data/julia_serialized/expression_deciles$(NEWDS ? "_new_ds" : "").json", "w") do file
    JSON.print(file, Dict([["$i" => decile for (i,decile) in enumerate(expression_deciles)]...; "singleton" => singleton_expr_vals]))
end

kw_test = KruskalWallisTest(expression_deciles...)

# Fig 3 statistics. Dunn's post-hoc tests between deciles are run by dunns_test.R.
suffix = NEWDS ? "_yn00" : ""
tests = DataFrame(Test = String[], Statistic = String[], Value = Float64[], N = Int[], p = Float64[])
push!(tests, ("Kruskal-Wallis, pair-average expression across dS deciles", "H", kw_test.H, sum(kw_test.n_i), pvalue(kw_test)))
push!(tests, ("Permutation test, Pearson correlation of pair-average expression with dS", "r", cor_test[2], nrow(paralog_data), cor_test[1]))

# Trend after the first decile: Spearman correlation of pair-average expression with dS over deciles 2-10
# (permutation test on ranks), and deciles 8-10 vs deciles 2-7
later = quantile_vals .>= 2
trend_test = perm_cor_2side(tiedrank(paralog_data.AvgExpr[later]), tiedrank(Vector(paralog_data.dS[later])))
push!(tests, ("Permutation test, Spearman correlation of pair-average expression with dS, deciles 2-10", "rho", trend_test[2], sum(later), trend_test[1]))
high_vs_mid = MannWhitneyUTest(paralog_data.AvgExpr[quantile_vals .>= 8], paralog_data.AvgExpr[2 .<= quantile_vals .<= 7])
push!(tests, ("Mann-Whitney U, pair-average expression, deciles 8-10 vs 2-7", "r", mwu_r(high_vs_mid), sum(quantile_vals .>= 2), pvalue(high_vs_mid)))
write_stats("expression_vs_ds$suffix", tests)

# All paralog pairs vs singletons
write_stats("expression_vs_ds_singletons$suffix",
    DataFrame([mwu_row("Pair-average expression of paralogs vs singleton expression", paralog_data.AvgExpr, singleton_expr_vals)]))

write_stats("expression_vs_ds_deciles$suffix", DataFrame(
    Decile = eachindex(expression_deciles),
    N = length.(expression_deciles),
    dS_min = [minimum(paralog_data.dS[quantile_vals .== q]) for q in eachindex(expression_deciles)],
    dS_max = [maximum(paralog_data.dS[quantile_vals .== q]) for q in eachindex(expression_deciles)],
    MedianExpr = median.(expression_deciles)))

# Plot expression vs dS (Figure 3)
fig_expr = plot([box_plt(paralog_data.AvgExpr[quantile_vals .== q], "$q") for q in sort(unique(quantile_vals), rev=false)], merge(Layout(title="Expression vs. 𝑑𝑆"), box_layout_expr))
savefig(fig_expr, "../../dicty_data/saved_figures/expression_vs_ds$suffix.html")
isinteractive() && display(fig_expr)
isinteractive() && display(plot(scatter(x=rollmean(paralog_data[:,"dS"], 100), y=rollmean(paralog_data.AvgExpr, 100), mode="markers", text=paralog_data.GeneID, marker_size=5), Layout(xaxis=attr(title="dS"), yaxis=attr(title="log2(expression)", range=[0, 4]))))

# Plot expression difference vs dS
isinteractive() && display(plot([box_plt(paralog_data.Diff[quantile_vals .== q], "$q") for q in sort(unique(quantile_vals), rev=false)], merge(Layout(title="Expression Difference"), box_layout_diff)))
isinteractive() && display(plot(scatter(x=rollmean(paralog_data[:,"dS"], 100), y=rollmean(paralog_data.Diff, 100), mode="markers", text=paralog_data.GeneID, marker_size=5), Layout(xaxis=attr(title="dS"), yaxis=attr(title="Expression difference"))))
CSV.write("../../dicty_data/filtered/paralog_expr_ds$(NEWDS ? "_new_ds" : "").csv", paralog_data)

df_adj = findall(row -> row.AvgExprAdj != Inf, eachrow(paralog_data))
quantile_vals_adj = quantile_vals[df_adj]
df_adj = paralog_data[df_adj,:]
adj_deciles = [log2.(df_adj.AvgExprAdj[quantile_vals_adj .== q] .+ 0.5) for q in sort(unique(quantile_vals_adj))]
fig_adj = plot([box_plt(adj_deciles[q], "$q") for q in sort(unique(quantile_vals_adj))],
        merge(Layout(title="Adjusted Expression vs 𝑑𝑆"), box_layout_expr))
savefig(fig_adj, "../../dicty_data/saved_figures/expression_vs_ds_adjusted$suffix.html")
isinteractive() && display(fig_adj)
if !NEWDS  # Figure S10
    ds_range(q) = @sprintf("%.2f–%.2f", extrema(paralog_data.dS[quantile_vals .== q])...)
    save_figure(plot([fig_box(adj_deciles[i], ds_range(q), BLUE) for (i, q) in enumerate(sort(unique(quantile_vals_adj)))],
                     fig_layout(xaxis=fig_axis(title="<i>dS</i> decile", tickangle=45, showgrid=false),
                                yaxis=fig_axis(title="Adjusted log<sub>2</sub>(TPM + 0.5)"))), "figS10"; width=800, height=550)
end

isinteractive() && display(plot(scatter(x=rollmean(df_adj[:,"dS"], 100),
                     y=rollmean(log2.(df_adj.AvgExprAdj .+ 0.5), 100), 
                     mode="markers", 
                     text=df_adj.GeneID, 
                     marker_size=5), 
            Layout(xaxis=attr(title="dS"), yaxis=attr(title="log2(expression)", range=[0, 4]))))
