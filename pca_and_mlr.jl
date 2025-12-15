using CategoricalArrays
using CSV
using DataFrames
using GLM
using LinearAlgebra
using MultivariateStats
using MultipleTesting
using NaturalSort
using PlotlyJS
using StatsBase
using YeoJohnsonTrans

data_dir = "../../dicty_data/"
id_method = "dS" # "dS", "dN", or "id"
id_threshold = 3

# Custom lib src:
include("./custom_lib/load_gff.jl")
include("./custom_lib/genomic_data.jl")
include("./custom_lib/enrichment_utils.jl")
include("./custom_lib/misc_utils.jl")

# helper function
function normalize(u::Vector{Float64})
    v = YeoJohnsonTrans.transform(u)
    return zscore(v)
end

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
insertcols!(expr_data, :Avg => mean.(eachrow(expr_data[:, 2:end])))
select!(expr_data, ["GeneID", "Avg"])

# log-transform data with an added pseudocount of 0.5
expr_data.Avg = normalize(expr_data.Avg)

# filter expression data to only include genes in the final gene list
filter!(row -> row.GeneID in final_gene_list, expr_data)

# Load peak data
peak_files = [filter(fn -> !contains(fn, r"_S[AB]_"), readdir(chip_peak_file_dir, join=true)); readdir(atac_peak_file_dir, join=true)]
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
for i in 1:nrow(full_df)
    gid = full_df.GeneID[i]
    pid = full_df.ParalogID[i]
    
        # H3K27ac
    gid_sig_h3k27ac = mean([mean(getsiginrange(get(ref_genome, gid), GeneRange(TSS(), TES(), 0, 0), ind)) for ind in [1, 4, 7]])
    pid_sig_h3k27ac = mean([mean(getsiginrange(get(ref_genome, pid), GeneRange(TSS(), TES(), 0, 0), ind)) for ind in [1, 4, 7]])
    full_df[i, :H3K27ac] = mean([gid_sig_h3k27ac, pid_sig_h3k27ac])

        # H3K4me3
    gid_sig_h3k4me3 = mean([mean(getsiginrange(get(ref_genome, gid), GeneRange(TSS(), TES(), 0, 0), ind)) for ind in [2, 5, 8]])
    pid_sig_h3k4me3 = mean([mean(getsiginrange(get(ref_genome, pid), GeneRange(TSS(), TES(), 0, 0), ind)) for ind in [2, 5, 8]])
    full_df[i, :H3K4me3] = mean([gid_sig_h3k4me3, pid_sig_h3k4me3])

        # H3K9me3
    gid_sig_h3k9me3 = mean([mean(getsiginrange(get(ref_genome, gid), GeneRange(TSS(), TES(), 0, 0), ind)) for ind in [3, 6, 9]])
    pid_sig_h3k9me3 = mean([mean(getsiginrange(get(ref_genome, pid), GeneRange(TSS(), TES(), 0, 0), ind)) for ind in [3, 6, 9]])
    full_df[i, :H3K9me3] = mean([gid_sig_h3k9me3, pid_sig_h3k9me3])

        # ATAC
    gid_sig_atac = mean([mean(getsiginrange(get(ref_genome, gid), GeneRange(TSS(), TES(), 0, 0), ind)) for ind in [10, 11, 12]])
    pid_sig_atac = mean([mean(getsiginrange(get(ref_genome, pid), GeneRange(TSS(), TES(), 0, 0), ind)) for ind in [10, 11, 12]])
    full_df[i, :ATAC] = mean([gid_sig_atac, pid_sig_atac])
end

display(plot([box(y= full_df.AvgExpr[full_df.H3K9me3 .> 0], name="w/ H3K9me3 (n = $(sum(full_df.H3K9me3 .> 0)))"), box(y= full_df.AvgExpr[full_df.H3K9me3 .<= 0], name="w/o H3K9me3 (n = $(sum(full_df.H3K9me3 .<= 0)))")], Layout(title="AvgExpr (w/ H3K9me3 vs w/o H3K9me3)")))
display(plot([box(y= full_df.AvgExpr[full_df.H3K27ac .> 0], name="w/ H3K27ac (n = $(sum(full_df.H3K27ac .> 0)))"), box(y= full_df.AvgExpr[full_df.H3K27ac .<= 0], name="w/o H3K27ac (n = $(sum(full_df.H3K27ac .<= 0)))")], Layout(title="AvgExpr (w/ H3K27ac vs w/o H3K27ac)")))
display(plot([box(y= full_df.AvgExpr[full_df.H3K4me3 .> 0], name="w/ H3K4me3 (n = $(sum(full_df.H3K4me3 .> 0)))"), box(y= full_df.AvgExpr[full_df.H3K4me3 .<= 0], name="w/o H3K4me3 (n = $(sum(full_df.H3K4me3 .<= 0)))")], Layout(title="AvgExpr (w/ H3K4me3 vs w/o H3K4me3)")))
display(plot([box(y= full_df.AvgExpr[full_df.ATAC .> 0], name="w/ ATAC (n = $(sum(full_df.ATAC .> 0)))"), box(y= full_df.AvgExpr[full_df.ATAC .<= 0], name="w/o ATAC (n = $(sum(full_df.ATAC .<= 0)))")], Layout(title="AvgExpr (w/ ATAC vs w/o ATAC)")))

# PCA ANALYSIS: #
full_df_backup = copy(full_df)  # Backup before z-scoring

# Z-score normalize the data
full_df = copy(full_df_backup)
full_df.AvgExpr = zscore(full_df.AvgExpr) # already log-transformed
full_df.H3K27ac = normalize(full_df.H3K27ac)
full_df.H3K4me3 = normalize(full_df.H3K4me3)
full_df.H3K9me3 = normalize(full_df.H3K9me3)
full_df.ATAC = normalize(full_df.ATAC)
full_df.dS = normalize(full_df.dS)

#--------------- PCA ---------------#
input_data = Matrix(full_df[:,3:end-1])'
pca_model = fit(PCA, input_data, maxoutdim=2)
pca_data = MultivariateStats.predict(pca_model, input_data)

# Use the loadings from the PCA to draw a biplot
loadings_matrix = MultivariateStats.loadings(pca_model)  # Variables x PCs
variable_names = ["AvgExpr", "H3K27ac", "H3K4me3", "H3K9me3", "ATAC"]

# Scale factor to make loadings visible relative to data spread
scale_factor = 2.0

# Create biplot traces
biplot_traces = GenericTrace[]

# Scatter plot of PCA scores
push!(biplot_traces, scatter(
    x=pca_data[1, :], y=pca_data[2, :],
    mode="markers",
    marker=attr(color=full_df.dS, opacity=0.5, size=5, colorscale="Viridis", colorbar=attr(title="dS")),
    name="Data points",
    showlegend=false
))

# Colors for each loading vector
colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]

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
    title="PCA Biplot",
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
savefig(biplot_fig, joinpath(data_dir, "pca_biplot.html"))

#--------------- MLR ---------------#
# For the ChIP/ATAC peaks, create binary variables indicating presence/absence of a peak
# and normalize the enrichment values for only those genes that have peaks
full_df = copy(full_df_backup)
insertcols!(full_df, :HasH3K27ac => full_df.H3K27ac .> 0)
insertcols!(full_df, :HasH3K4me3 => full_df.H3K4me3 .> 0)
insertcols!(full_df, :HasH3K9me3 => full_df.H3K9me3 .> 0)
insertcols!(full_df, :HasATAC => full_df.ATAC .> 0)
full_df.H3K27ac[full_df.H3K27ac .> 0] = normalize(full_df.H3K27ac[full_df.H3K27ac .> 0])
full_df.H3K4me3[full_df.H3K4me3 .> 0] = normalize(full_df.H3K4me3[full_df.H3K4me3 .> 0])
full_df.H3K9me3[full_df.H3K9me3 .> 0] = normalize(full_df.H3K9me3[full_df.H3K9me3 .> 0])
full_df.ATAC[full_df.ATAC .> 0] = normalize(full_df.ATAC[full_df.ATAC .> 0])
full_df.AvgExpr = normalize(full_df.AvgExpr)
full_df.dS = normalize(full_df.dS)

# calculate the VIFs by calculating the correlation matrix and taking the diagonal of its inverse
mat = Matrix(full_df[:,3:end])
cor_mat = cor(mat)
inv_cor_mat = inv(cor_mat)
vif_values = diag(inv_cor_mat)

mlr_model = lm(
    @formula(
        dS ~ 
        HasH3K27ac + H3K27ac + 
        HasH3K4me3 + H3K4me3 + 
        HasH3K9me3 + H3K9me3 + 
        HasATAC + ATAC + 
        AvgExpr
    ), 
    full_df)
mlr_table = DataFrame(coeftable(mlr_model))
rename!(mlr_table, [:Predictor, :Coef, :StdErr, :tvalue, :pvalue, :Lower95CI, :Upper95CI])
r_squared = adjr²(mlr_model)
mlr_table.pvalue_adj = adjust(mlr_table[!, :pvalue], BenjaminiHochberg())

sort!(mlr_table, :pvalue_adj)
CSV.write(joinpath(data_dir, "mlr_results.tsv"), mlr_table)