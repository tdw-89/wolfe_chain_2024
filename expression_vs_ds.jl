include("prelude.jl")

using CategoricalArrays
using NaturalSort
using RollingFunctions
using JSON


using .EnrichmentUtils
using .MiscUtils

function avg_enrich(pair, inds, gene_range)
    range_1 = deepcopy(gene_range)
    range_2 = deepcopy(gene_range)
    avgs = []

    for ind in inds
        got_sig_1 = false
        got_sig_2 = false
        sig_1 = nothing
        sig_2 = nothing

        while !got_sig_1
            
            sig_1 = getsiginrange(pair[1], range_1, ind)
            if ismissing(sig_1)
                range_1 = GeneRange(TSS(), TES(), 0, 0)
            else
                got_sig_1 = true
            end
        end

        while !got_sig_2
            sig_2 = getsiginrange(pair[2], range_2, ind)
            if ismissing(sig_2)
                range_2 = GeneRange(TSS(), TES(), 0, 0)
            else
                got_sig_2 = true
            end
        end
        ind_avg = mean([mean(sig_1), mean(sig_2)])
        push!(avgs, ind_avg)
    end

    return mean(avgs)
end

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

# # Load final filtered gene list
# final_gene_list = open(final_gene_list) do file
#     readlines(file)
# end

# Load expression data
expr_data = CSV.read(expr_data_file, DataFrame)

# Rename the expression columns to a letter indicating with life-cycle stage the sample
# was taken from.
rename!(expr_data, ["GeneID", "V", "S", "M", "F"])
expr_data[!,2:end] = log2.(expr_data[!,2:end] .+ 0.5)

# Load peak data
peak_files = [filter(fn -> !contains(fn, r"_S[AB]_"), readdir(chip_peak_file_dir, join=true)); readdir(atac_peak_file_dir, join=true)]
filter!(fn -> endswith(fn, ".narrowPeak"), peak_files)
peak_data = binpeaks(peak_files, chrom_lengths_file)

# Load paralog and singleton data
paralog_data = CSV.read(paralog_file, DataFrame)
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
singleton_expr_vals = [gene.rnas[1].expression[1] for gene in singletons if !ismissing(gene) && !isempty(gene.rnas)]
expression_deciles = [paralog_data.AvgExpr[quantile_vals .== q] for q in sort(unique(quantile_vals))]

open("../../dicty_data/julia_serialized/expression_deciles.json", "w") do file
    JSON.print(file, Dict([["$i" => decile for (i,decile) in enumerate(expression_deciles)]...; "singleton" => singleton_expr_vals]))
end

kw_test = KruskalWallisTest(expression_deciles...)
mwu_tests = [MannWhitneyUTest(expression_decile, singleton_expr_vals) for expression_decile in expression_deciles]
adj_pvals = adjust(pvalue.(mwu_tests), BenjaminiHochberg())

# Get the significance of the overall difference
dup_vs_single = pvalue(MannWhitneyUTest(singleton_expr_vals, 
                                        paralog_data.AvgExpr))

# Get correlation between pair-average expression and dS
expr_ds_corr = cor(paralog_data.AvgExpr, paralog_data[!,3]), pvalue(CorrelationTest(paralog_data.AvgExpr, paralog_data[!,3]))

display(plot([box_plt(paralog_data.AvgExpr[quantile_vals .== q], "$q") for q in sort(unique(quantile_vals), rev=false)], 
                merge(Layout(title="Expression vs. 𝑑𝑆"), box_layout_expr)))
display(plot(scatter(x=rollmean(paralog_data[:,"dS"], 100), y=rollmean(paralog_data.AvgExpr, 100), mode="markers", text=paralog_data.GeneID, marker_size=5), Layout(xaxis=attr(title="dS"), yaxis=attr(title="log2(expression)", range=[0, 4]))))

display(plot([box_plt(paralog_data.Diff[quantile_vals .== q], "$q") for q in sort(unique(quantile_vals), rev=false)], 
                merge(Layout(title="Expression Difference"), box_layout_diff)))
display(plot(scatter(x=rollmean(paralog_data[:,"dS"], 100), y=rollmean(paralog_data.Diff, 100), mode="markers", text=paralog_data.GeneID, marker_size=5), Layout(xaxis=attr(title="dS"), yaxis=attr(title="Expression difference"))))

df_adj = findall(row -> row.AvgExprAdj != Inf, eachrow(paralog_data))
quantile_vals_adj = quantile_vals[df_adj]
df_adj = paralog_data[df_adj,:]
quantile_text = parse_quantile(sort(unique(String.(quantile_labels)), lt=natural), digs=2)
adj_deciles = [log.(df_adj.AvgExprAdj[quantile_vals_adj .== q] .+ 0.5) for q in sort(unique(quantile_vals_adj))]
display(plot([box_plt(adj_deciles[q], "$q") for q in sort(unique(quantile_vals_adj))], 
        merge(Layout(title="Adjusted Expression vs 𝑑𝑆"), box_layout_expr)))

display(plot(scatter(x=rollmean(df_adj[:,"dS"], 100),
                     y=rollmean(log.(df_adj.AvgExprAdj .+ 0.5), 100), 
                     mode="markers", 
                     text=df_adj.GeneID, 
                     marker_size=5), 
            Layout(xaxis=attr(title="dS"), yaxis=attr(title="log2(expression)", range=[0, 4]))))