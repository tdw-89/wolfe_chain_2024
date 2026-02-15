include("prelude.jl")

using .EnrichmentUtils
import .MiscUtils: normalize_yj

using CategoricalArrays
using GLM
using LinearAlgebra
using MultivariateStats
using NaturalSort

const ϵ = 0.001

# Helper functions
function get_avg_signal(gene::Gene, idx::Int64)
    sig = getsiginrange(gene, GeneRange(TSS(), TES(), -500, 500), idx)
    return mean(sig)
end

function distance(geneA::Gene, geneB::Gene)
    ismissing(geneA.scaffold) && error("Gene $(geneA.id) is missing scaffold information")
    ismissing(geneB.scaffold) && error("Gene $(geneB.id) is missing scaffold information")
    geneA.scaffold.name != geneB.scaffold.name && return Inf
    if geneA.gene_end < geneB.gene_start
        return geneB.gene_start - geneA.gene_end
    elseif geneB.gene_end < geneA.gene_start
        return geneA.gene_start - geneB.gene_end
    else
        return 0
    end
end

function logit(vec::Vector{Float64})
    vec = clamp.(vec, ϵ, 1 - ϵ)
    return log.(vec ./ (1 .- vec))
end

function logit_z(vec::Vector{Float64})
    log_odd_vec = logit(vec)
    return zscore(log_odd_vec)
end

data_dir = "../../dicty_data/"
life_cycle = "M" # "F" "M" or "All"
PCA_analysis = true

k27ac_upstream, k4me3_upstream, k9me3_upstream, atac_upstream = -500, -500, -500, -500
k27ac_downstream, k4me3_downstream, k9me3_downstream, atac_downstream = 500, 500, 500, 500

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

# Rename the expression columns to a letter indicating which life-cycle stage the sample
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
expr_data.Avg = log2.(expr_data.Avg .+ 0.5)

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
filtered_gene_list = [gene for gene in ref_genome.genes[2] if gene.id in final_gene_list]

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

idx_groups = [idxs_h3k27ac, idxs_h3k4me3, idxs_h3k9me3, idxs_atac]
group_names = ["H3K27ac", "H3K4me3", "H3K9me3", "ATAC"]
genes_means = [Float64[] for _ in 1:length(idx_groups)]
group_means = zeros(length(idx_groups))

for (i, idx_group) in enumerate(idx_groups)
    for gene in filtered_gene_list
        for idx in idx_group
            x = get_avg_signal(gene, idx)
            push!(genes_means[i], x)
        end
    end
    group_means[i] = mean(genes_means[i])
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

# Calculate the distance between paralogs
paralog_dist = Float64[]
for i in 1:nrow(full_df)
    gid = full_df.GeneID[i]
    pid = full_df.ParalogID[i]
    gene1 = get(ref_genome, gid)
    gene2 = get(ref_genome, pid)
    if gene1.scaffold.name == gene2.scaffold.name
        push!(paralog_dist, distance(gene1, gene2))
    else
        push!(paralog_dist, Inf)
    end
end
full_df.ParalogDist = paralog_dist

full_df_backup = copy(full_df) # Backup before normalization

# PCA ANALYSIS: #
full_df = copy(full_df_backup)
# Create indicator variables
full_df.HasH3K27ac = ifelse.(full_df.H3K27ac .> 0, 1, 0)
full_df.HasH3K4me3 = ifelse.(full_df.H3K4me3 .> 0, 1, 0)
full_df.HasH3K9me3 = ifelse.(full_df.H3K9me3 .> 0, 1, 0)
full_df.HasATAC = ifelse.(full_df.ATAC .> 0, 1, 0)
full_df.SameChrom = ifelse.(isfinite.(full_df.ParalogDist), 1, 0)
full_df.ParalogDist[isinf.(full_df.ParalogDist)] .= 0

# Normalize the data
full_df.H3K27ac[full_df.H3K27ac .> 0] = logit_z(full_df.H3K27ac[full_df.H3K27ac .> 0])
full_df.HasH3K27ac = zscore(full_df.HasH3K27ac)
full_df.H3K4me3[full_df.H3K4me3 .> 0] = logit_z(full_df.H3K4me3[full_df.H3K4me3 .> 0])
full_df.HasH3K4me3 = zscore(full_df.HasH3K4me3)
full_df.H3K9me3[full_df.H3K9me3 .> 0] = logit_z(full_df.H3K9me3[full_df.H3K9me3 .> 0])
full_df.HasH3K9me3 = zscore(full_df.HasH3K9me3)
full_df.ATAC[full_df.ATAC .> 0] = logit_z(full_df.ATAC[full_df.ATAC .> 0])
full_df.HasATAC = zscore(full_df.HasATAC)
full_df.ParalogDist[full_df.ParalogDist .> 0] = log.(full_df.ParalogDist[full_df.ParalogDist .> 0]) |> zscore
full_df.SameChrom = zscore(full_df.SameChrom)
full_df.AvgExpr = zscore(full_df.AvgExpr)
full_df.dS = log.(full_df.dS .+ ϵ) |> zscore

#--------------- PCA ---------------#
input_data = Matrix(full_df[:,4:end])'
pca_model = fit(PCA, input_data, maxoutdim=size(input_data, 1))
pca_data = MultivariateStats.predict(pca_model, input_data)

# Use the loadings from the PCA to draw a biplot
loadings_matrix = MultivariateStats.loadings(pca_model) # Variables x PCs
variable_names = names(full_df)[4:end]

# Scale factor to make loadings visible relative to data spread
scale_factor = 3.0

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

# Add loading vectors as lines with arrowhead markers
for (i, var_name) in enumerate(variable_names)
    pc1_loading = loadings_matrix[i, 1] * scale_factor
    pc2_loading = loadings_matrix[i, 2] * scale_factor
    
    # Add the vector line (from origin to loading)
    push!(biplot_traces, scatter(
        x=[0, pc1_loading],
        y=[0, pc2_loading],
        mode="lines+markers",
        line=attr(width=2),
        marker=attr(symbol="arrow", size=12, angleref="previous"),
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
savefig(biplot_fig, joinpath(data_dir, "pca_biplot_$(life_cycle).html"))

#--------------- MLR on PCA components ---------------#
pca_df = DataFrame(PC1 = pca_data[1, :],
                    PC2 = pca_data[2, :],
                    PC3 = pca_data[3, :],
                    PC4 = pca_data[4, :],
                    PC5 = pca_data[5, :],
                    PC6 = pca_data[6, :],
                    PC7 = pca_data[7, :],
                    PC8 = pca_data[8, :],
                    PC9 = pca_data[9, :],
                    PC10 = pca_data[10, :],
                    dS = full_df.dS)

mlr_pca_model = lm(@formula(dS ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10), pca_df)
mlr_pca_table = DataFrame(coeftable(mlr_pca_model))
r2_squared_pca = r²(mlr_pca_model)
rename!(mlr_pca_table, [:Predictor, :Coef, :StdErr, :tvalue, :pvalue, :Lower95CI, :Upper95CI])
mlr_pca_table.pvalue_adj = adjust(mlr_pca_table[!, :pvalue], BenjaminiHochberg())
mlr_pca_table = mlr_pca_table[2:end, :]  # Remove intercept row (p value of 1)
mlr_pca_table.prop_var_explained = pca_model.prinvars ./ pca_model.tprinvar

# For each of the original input variables, add a column with the loading on each PC
for (i, var_name) in enumerate(variable_names)
    mlr_pca_table[!, Symbol(var_name * "_loading")] = loadings_matrix[i, :]
end

CSV.write(joinpath(data_dir, "mlr_pca_results_$(life_cycle).csv"), mlr_pca_table)

# Plot results:
plot_table = select(mlr_pca_table, [1, 2, 8:ncol(mlr_pca_table)...])
DataFrames.transform!(plot_table, names(plot_table)[5:end] .=> (col -> col .* sign.(plot_table.Coef)) .=> names(plot_table)[5:end])
for i in 1:nrow(plot_table)
    loadings = Vector(plot_table[i,5:end])
    loadings_total = abs.(loadings) |> sum
    plot_table[i,5:end] .= abs.(loadings) ./ loadings_total .* plot_table.prop_var_explained[i] .* 100
end

#--------------- Stacked Bar Chart ---------------#
# Get loading column names (columns 5 to end)
loading_col_names = names(plot_table)[5:end]

# Color mapping: related variables (coverage + indicator) share the same color
color_map = Dict(
    "AvgExpr"      => "#1f77b4",  # blue
    "H3K27ac"      => "#ff7f0e",  # orange
    "HasH3K27ac"   => "#ff7f0e",  # orange
    "H3K4me3"      => "#2ca02c",  # green
    "HasH3K4me3"   => "#2ca02c",  # green
    "H3K9me3"      => "#d62728",  # red
    "HasH3K9me3"   => "#d62728",  # red
    "ATAC"         => "#9467bd",  # purple
    "HasATAC"      => "#9467bd",  # purple
    "ParalogDist"  => "#8c564b",  # brown
    "SameChrom"    => "#8c564b",  # brown
)

# Create bar traces for each loading variable
bar_traces = GenericTrace[]
for (j, col_name) in enumerate(loading_col_names)
    # Clean up the display name by removing "_loading" suffix
    display_name = replace(col_name, "_loading" => "")
    
    # Get color from map, fallback to gray if not found
    marker_color = get(color_map, display_name, "#7f7f7f")
    
    push!(bar_traces, bar(
        x = plot_table.Predictor,
        y = plot_table[!, col_name],
        name = display_name,
        marker = attr(color = marker_color),
        hovertemplate = "%{x}<br>$display_name: %{y:.2f}%<extra></extra>"
    ))
end

# Create layout with barmode="stack" for stacked bars
bar_layout = Layout(
    title = attr(
        text = "PCA Component Loading Contributions - $life_cycle<br><sub>R<sup>2</sup> = $(round(r2_squared_pca, digits=3))</sub>",
        x = 0.5,
        xanchor = "center"
    ),
    xaxis = attr(
        title = "Principal Component",
        tickmode = "array",
        tickvals = plot_table.Predictor,
        categoryorder = "array",
        categoryarray = plot_table.Predictor
    ),
    yaxis = attr(
        title = "Loading Proportion",
        gridcolor = "lightgray"
    ),
    barmode = "stack",  # Stacks all bars upward from zero
    legend = attr(
        x = 1.02,
        y = 1,
        xanchor = "left",
        yanchor = "top",
        title = attr(text = "Variables")
    ),
    margin = attr(r = 180, b = 80),
    plot_bgcolor = "white",
    paper_bgcolor = "white",
    hovermode = "closest"
)

bar_fig = plot(bar_traces, bar_layout)
display(bar_fig)
savefig(bar_fig, joinpath(data_dir, "pca_contributions_stacked_$(life_cycle).html"))