include("prelude.jl")

using .EnrichmentUtils
import .MiscUtils: normalize_yj

using CategoricalArrays
using GLM
using LinearAlgebra
using MultivariateStats
using NaturalSort

data_dir = "../../dicty_data/"
life_cycle = "All" # "F" "M" "V" or "All"
sigregions_only = false

k27ac_upstream, k4me3_upstream, k9me3_upstream, atac_upstream = -500, -500, -500, -500
k27ac_downstream, k4me3_downstream, k9me3_downstream, atac_downstream = 500, 500, 500, 500
if sigregions_only
    k27ac_upstream, k27ac_downstream = -140, 148
    k4me3_upstream, k4me3_downstream = -200, 164
    k9me3_upstream, k9me3_downstream = -500, 500
    atac_upstream, atac_downstream = -500, 219
end

id_threshold = 3

# Genome data
gff_data = joinpath(data_dir, "AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3")
chrom_lengths_file = joinpath(data_dir, "AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt")
final_gene_list = "../../dicty_data/filtered/final_gene_list.txt"

# Expression data
expr_data_file = joinpath(data_dir, "filtered/expr_data_filt_kallisto_ensembl52_single.tsv")

# Peak files
chip_peak_file_dir = joinpath(data_dir, "wang_et_al/processed/run_1_ensembl52/")
atac_peak_file_dir = joinpath(data_dir, "wang_et_al/processed/run_2_ensembl52/")

# Paralog data
paralog_file = joinpath(data_dir, "filtered/paralog_filt.tsv")

# Load final filtered gene list
final_gene_list = open(final_gene_list) do file
    readlines(file)
end

# Load expression data
expr_data = CSV.read(expr_data_file, DataFrame)

# Rename the expression columns to a letter indicating with life-cycle stage the sample
# was taken from.
rename!(expr_data, ["GeneID", "V", "S", "M", "F"])

# Remove 'S' stage expression data, and rearrange so that the samples are in the order they
# will be for the peak data
select!(expr_data, ["GeneID", "F", "M", "V"])

# Average gene expression counts across all life-cycle stages
if life_cycle != "All"
    select!(expr_data, ["GeneID", life_cycle])
end
insertcols!(expr_data, :Avg => mean.(eachrow(expr_data[:, 2:end])))
select!(expr_data, ["GeneID", "Avg"])

# log-transform data with an added pseudocount of 0.5
expr_data.Avg = normalize_yj(expr_data.Avg)

# filter expression data to only include genes in the final gene list
filter!(row -> row.GeneID in final_gene_list, expr_data)

# Load peak data
peak_files = [filter(fn -> !contains(fn, r"_S[AB]_"), readdir(chip_peak_file_dir, join=true)); readdir(atac_peak_file_dir, join=true)]
if life_cycle != "All"
    peak_files = filter(fn -> contains(fn, Regex("_$life_cycle[AB]?_")), peak_files)
end
filter!(fn -> endswith(fn, ".narrowPeak"), peak_files)
peak_data = binpeaks(peak_files, chrom_lengths_file)

# Load paralog
paralog_data = CSV.read(paralog_file, DataFrame)

# Filter the paralog data to just those that have expression data
filter!(row -> row.GeneID in expr_data.GeneID && row.ParalogID in expr_data.GeneID, paralog_data)

# Filter paralogs according to dS threshold
select!(paralog_data, ["GeneID", "ParalogID", "dS"])
filter!(row -> row["dS"] <= id_threshold, paralog_data)

# Load reference genome object
ref_genome = loadgenome(gff_data, chrom_lengths_file)

# Add data to the reference
addexpression!(ref_genome, expr_data)
addtogenes!(ref_genome, peak_data)

# Create a dataframe for all paralog data
full_df = copy(paralog_data)

# Get the average expression for each gene
insertcols!(full_df, :AvgExpr => zeros(nrow(full_df)))
for (i, (gid, pid)) in enumerate(zip(full_df.GeneID, full_df.ParalogID))
    gid_expr = only(getexpression(get(ref_genome, gid)))
    pid_expr = only(getexpression(get(ref_genome, pid)))
    full_df[i, :AvgExpr] = mean([gid_expr, pid_expr])
end

# Get the average enrichment for each gene for each mark
insertcols!(full_df, :H3K27ac => zeros(nrow(full_df)))
insertcols!(full_df, :H3K4me3 => zeros(nrow(full_df)))
insertcols!(full_df, :H3K9me3 => zeros(nrow(full_df)))
insertcols!(full_df, :ATAC => zeros(nrow(full_df)))

# indexes for each mark in each tissue
idxs_h3k27ac = [1, 4, 7]
idxs_h3k4me3 = [2, 5, 8]
idxs_h3k9me3 = [3, 6, 9]
idxs_atac = [10, 11, 12]
if life_cycle != "All"
    idxs_h3k27ac = [1]
    idxs_h3k4me3 = [2]
    idxs_h3k9me3 = [3]
    idxs_atac = [4]
end
for i in 1:nrow(full_df)
    gid = full_df.GeneID[i]
    pid = full_df.ParalogID[i]
    
        # H3K27ac
    gid_sig_h3k27ac = mean([mean(getsiginrange(get(ref_genome, gid), GeneRange(TSS(), TES(), k27ac_upstream, k27ac_downstream), ind)) for ind in idxs_h3k27ac])
    pid_sig_h3k27ac = mean([mean(getsiginrange(get(ref_genome, pid), GeneRange(TSS(), TES(), k27ac_upstream, k27ac_downstream), ind)) for ind in idxs_h3k27ac])
    full_df[i, :H3K27ac] = mean([gid_sig_h3k27ac, pid_sig_h3k27ac])

        # H3K4me3
    gid_sig_h3k4me3 = mean([mean(getsiginrange(get(ref_genome, gid), GeneRange(TSS(), TES(), k4me3_upstream, k4me3_downstream), ind)) for ind in idxs_h3k4me3])
    pid_sig_h3k4me3 = mean([mean(getsiginrange(get(ref_genome, pid), GeneRange(TSS(), TES(), k4me3_upstream, k4me3_downstream), ind)) for ind in idxs_h3k4me3])
    full_df[i, :H3K4me3] = mean([gid_sig_h3k4me3, pid_sig_h3k4me3])

        # H3K9me3
    gid_sig_h3k9me3 = mean([mean(getsiginrange(get(ref_genome, gid), GeneRange(TSS(), TES(), k9me3_upstream, k9me3_downstream), ind)) for ind in idxs_h3k9me3])
    pid_sig_h3k9me3 = mean([mean(getsiginrange(get(ref_genome, pid), GeneRange(TSS(), TES(), k9me3_upstream, k9me3_downstream), ind)) for ind in idxs_h3k9me3])
    full_df[i, :H3K9me3] = mean([gid_sig_h3k9me3, pid_sig_h3k9me3])

        # ATAC
    gid_sig_atac = mean([mean(getsiginrange(get(ref_genome, gid), GeneRange(TSS(), TES(), atac_upstream, atac_downstream), ind)) for ind in idxs_atac])
    pid_sig_atac = mean([mean(getsiginrange(get(ref_genome, pid), GeneRange(TSS(), TES(), atac_upstream, atac_downstream), ind)) for ind in idxs_atac])
    full_df[i, :ATAC] = mean([gid_sig_atac, pid_sig_atac])
end

# display(plot([box(y= full_df.AvgExpr[full_df.H3K9me3 .> 0], name="w/ H3K9me3 (n = $(sum(full_df.H3K9me3 .> 0)))"), box(y= full_df.AvgExpr[full_df.H3K9me3 .<= 0], name="w/o H3K9me3 (n = $(sum(full_df.H3K9me3 .<= 0)))")], Layout(title="AvgExpr (w/ H3K9me3 vs w/o H3K9me3)")))
# display(plot([box(y= full_df.AvgExpr[full_df.H3K27ac .> 0], name="w/ H3K27ac (n = $(sum(full_df.H3K27ac .> 0)))"), box(y= full_df.AvgExpr[full_df.H3K27ac .<= 0], name="w/o H3K27ac (n = $(sum(full_df.H3K27ac .<= 0)))")], Layout(title="AvgExpr (w/ H3K27ac vs w/o H3K27ac)")))
# display(plot([box(y= full_df.AvgExpr[full_df.H3K4me3 .> 0], name="w/ H3K4me3 (n = $(sum(full_df.H3K4me3 .> 0)))"), box(y= full_df.AvgExpr[full_df.H3K4me3 .<= 0], name="w/o H3K4me3 (n = $(sum(full_df.H3K4me3 .<= 0)))")], Layout(title="AvgExpr (w/ H3K4me3 vs w/o H3K4me3)")))
# display(plot([box(y= full_df.AvgExpr[full_df.ATAC .> 0], name="w/ ATAC (n = $(sum(full_df.ATAC .> 0)))"), box(y= full_df.AvgExpr[full_df.ATAC .<= 0], name="w/o ATAC (n = $(sum(full_df.ATAC .<= 0)))")], Layout(title="AvgExpr (w/ ATAC vs w/o ATAC)")))

# PCA ANALYSIS: #
full_df_backup = copy(full_df)  # Backup before z-scoring

# Z-score normalize the data
full_df = copy(full_df_backup)
full_df.H3K27ac = normalize_yj(full_df.H3K27ac)
full_df.H3K4me3 = normalize_yj(full_df.H3K4me3)
full_df.H3K9me3 = normalize_yj(full_df.H3K9me3)
full_df.ATAC = normalize_yj(full_df.ATAC)
full_df.dS = normalize_yj(full_df.dS)

#--------------- PCA ---------------#
input_data = Matrix(full_df[:,3:end])'
pca_model = fit(PCA, input_data, maxoutdim=2)
pca_data = MultivariateStats.predict(pca_model, input_data)

# Use the loadings from the PCA to draw a biplot
loadings_matrix = MultivariateStats.loadings(pca_model) # Variables x PCs
variable_names = ["dS", "AvgExpr", "H3K27ac", "H3K4me3", "H3K9me3", "ATAC"]

# Scale factor to make loadings visible relative to data spread
scale_factor = 2.0

# Create biplot traces
biplot_traces = GenericTrace[]

# Scatter plot of PCA scores
push!(biplot_traces, scatter(
    x=pca_data[1, :], y=pca_data[2, :],
    mode="markers",
    marker=attr(size=6, color="rgba(0,0,0,0.5)"),
    name="Data points",
    showlegend=false
))

# Colors for each loading vector
colors = ["#8c564b", "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

# Add loading vectors as lines with arrowhead markers
for (i, var_name) in enumerate(variable_names)
    pc1_loading = loadings_matrix[i, 1] * scale_factor
    pc2_loading = loadings_matrix[i, 2] * scale_factor
    
    # Add the vector line (from origin to loading)
    push!(biplot_traces, scatter(
        x=[0, pc1_loading],
        y=[0, pc2_loading],
        mode="lines+markers",
        line=attr(color=colors[i], width=2),
        marker=attr(symbol="arrow", size=12, angleref="previous", color=colors[i]),
        name=var_name
    ))
end

# Create biplot layout
biplot_layout = Layout(
    title="PCA Biplot $life_cycle",
    xaxis_title="PC1",
    yaxis_title="PC2",
    showlegend=true,
    legend=attr(x=0, y=-0.15, xanchor="left", orientation="h"),
    xaxis=attr(zeroline=true, zerolinecolor="black", zerolinewidth=1),
    yaxis=attr(zeroline=true, zerolinecolor="black", zerolinewidth=1, scaleanchor="x", scaleratio=1),
    margin=attr(b=100)
)

biplot_fig = plot(biplot_traces, biplot_layout)
display(biplot_fig)
savefig(biplot_fig, joinpath(data_dir, "pca_biplot_$(life_cycle)$(sigregions_only ? "sigregions" : "allregions").html"))

#--------------- MLR ---------------#
# For the ChIP/ATAC peaks, create binary variables indicating presence/absence of a peak
# and normalize the enrichment values for only those genes that have peaks
full_df_mlr = copy(full_df_backup)
full_df_mlr.H3K27ac = normalize_yj(full_df_mlr.H3K27ac)
full_df_mlr.H3K4me3 = normalize_yj(full_df_mlr.H3K4me3)
full_df_mlr.H3K9me3 = normalize_yj(full_df_mlr.H3K9me3)
full_df_mlr.ATAC = normalize_yj(full_df_mlr.ATAC)
full_df_mlr.dS = normalize_yj(full_df_mlr.dS)

# calculate the VIFs by calculating the correlation matrix and taking the diagonal of its inverse
mat = Matrix(full_df_mlr[:,3:end])
cor_mat = cor(mat)
# Create a heatmap of the correlation matrix, with values annotated
heatmap_trace = heatmap(
    z=cor_mat,
    x=names(full_df_mlr)[3:end],
    y=names(full_df_mlr)[3:end],
    colorscale="Viridis",
    zmin=-1,
    zmax=1,
    colorbar=attr(title="Correlation")
)

heatmap_layout = Layout(
    title="Correlation Matrix Heatmap $life_cycle",
    xaxis=attr(tickangle=-45),
    yaxis=attr(autorange="reversed")
)

heatmap_fig = plot(heatmap_trace, heatmap_layout)
display(heatmap_fig)
savefig(heatmap_fig, joinpath(data_dir, "correlation_heatmap_$(life_cycle)$(sigregions_only ? "sigregions" : "allregions").html"))

inv_cor_mat = inv(cor_mat)
vif_values = diag(inv_cor_mat)

mlr_model = lm(
    @formula(
        dS ~
        H3K27ac +
        H3K4me3 +
        H3K9me3 +
        ATAC +
        AvgExpr
    ), 
    full_df_mlr)
mlr_table = DataFrame(coeftable(mlr_model))
rename!(mlr_table, [:Predictor, :Coef, :StdErr, :tvalue, :pvalue, :Lower95CI, :Upper95CI])
r_squared = r²(mlr_model)
mlr_table.pvalue_adj = adjust(mlr_table[!, :pvalue], BenjaminiHochberg())

sort!(mlr_table, :pvalue_adj)
CSV.write(joinpath(data_dir, "mlr_results_$(life_cycle)$(sigregions_only ? "sigregions" : "allregions").csv"), mlr_table)
stats_file = open(joinpath(data_dir, "mlr_stats_$(life_cycle)$(sigregions_only ? "sigregions" : "allregions").txt"), "w")
write(stats_file, "r²: $r_squared\n\nVIF Values:\n")
for (i, col_name) in enumerate(names(full_df_mlr)[3:end])
    write(stats_file, "$col_name: $(vif_values[i])\n")
end
close(stats_file)