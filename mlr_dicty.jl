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
function normalize_expression(df::DataFrame)
    num_obs = nrow(df)
    expr_cols = names(df)[2:end]
    pool_expr = reduce(vcat, [df[!, col] for col in expr_cols])
    pool_norm = normalize_yj(pool_expr)
    expr_mat = reshape(pool_norm, num_obs, length(expr_cols))
    new_df = DataFrame(expr_mat, Symbol.(expr_cols))
    insertcols!(new_df, 1, :GeneID => df.GeneID)
    return new_df
end

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
life_cycle = "All" # "F" "M" "V" or "All"

gene_range = GeneRange(TSS(), TES(), -500, 500)

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

# Average gene expression counts across all life-cycle stages
if life_cycle != "All"
    select!(expr_data, ["GeneID", life_cycle])
    expr_data[!, :Avg] = expr_data[!, life_cycle]
    select!(expr_data, ["GeneID", "Avg"])
    expr_data.Avg = normalize_yj(expr_data.Avg)
else
    expr_data = normalize_expression(expr_data)
    insertcols!(expr_data, :Avg => mean.(eachrow(expr_data[:, 2:end])))
    select!(expr_data, ["GeneID", "Avg"])
end

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

for i in 1:nrow(full_df)
    gid = full_df.GeneID[i]
    pid = full_df.ParalogID[i]
    
        # H3K27ac
    gid_sig_h3k27ac = mean([mean(getsiginrange(get(ref_genome, gid), gene_range, ind)) for ind in idxs_h3k27ac])
    pid_sig_h3k27ac = mean([mean(getsiginrange(get(ref_genome, pid), gene_range, ind)) for ind in idxs_h3k27ac])
    full_df[i, :H3K27ac] = mean([gid_sig_h3k27ac, pid_sig_h3k27ac])

        # H3K4me3
    gid_sig_h3k4me3 = mean([mean(getsiginrange(get(ref_genome, gid), gene_range, ind)) for ind in idxs_h3k4me3])
    pid_sig_h3k4me3 = mean([mean(getsiginrange(get(ref_genome, pid), gene_range, ind)) for ind in idxs_h3k4me3])
    full_df[i, :H3K4me3] = mean([gid_sig_h3k4me3, pid_sig_h3k4me3])

        # H3K9me3
    gid_sig_h3k9me3 = mean([mean(getsiginrange(get(ref_genome, gid), gene_range, ind)) for ind in idxs_h3k9me3])
    pid_sig_h3k9me3 = mean([mean(getsiginrange(get(ref_genome, pid), gene_range, ind)) for ind in idxs_h3k9me3])
    full_df[i, :H3K9me3] = mean([gid_sig_h3k9me3, pid_sig_h3k9me3])

        # ATAC
    if life_cycle != "V"        
        gid_sig_atac = mean([mean(getsiginrange(get(ref_genome, gid), gene_range, ind)) for ind in idxs_atac])
        pid_sig_atac = mean([mean(getsiginrange(get(ref_genome, pid), gene_range, ind)) for ind in idxs_atac])
        full_df[i, :ATAC] = mean([gid_sig_atac, pid_sig_atac])
    end
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

#--------------- MLR ---------------#
# Create indicator variables
full_df.HasH3K27ac = ifelse.(full_df.H3K27ac .> 0, 1, 0)
full_df.HasH3K4me3 = ifelse.(full_df.H3K4me3 .> 0, 1, 0)
full_df.HasH3K9me3 = ifelse.(full_df.H3K9me3 .> 0, 1, 0)
full_df.HasATAC = ifelse.(full_df.ATAC .> 0, 1, 0)
full_df.SameChrom = ifelse.(isfinite.(full_df.ParalogDist), 1, 0)

# Normalize the data
full_df.H3K27ac[full_df.H3K27ac .> 0] = logit_z(full_df.H3K27ac[full_df.H3K27ac .> 0])
full_df.H3K4me3[full_df.H3K4me3 .> 0] = logit_z(full_df.H3K4me3[full_df.H3K4me3 .> 0])
full_df.H3K9me3[full_df.H3K9me3 .> 0] = logit_z(full_df.H3K9me3[full_df.H3K9me3 .> 0])
full_df.ATAC[full_df.ATAC .> 0] = logit_z(full_df.ATAC[full_df.ATAC .> 0])
full_df.ParalogDist[isfinite.(full_df.ParalogDist)] = normalize_yj(full_df.ParalogDist[isfinite.(full_df.ParalogDist)])
full_df.ParalogDist[isinf.(full_df.ParalogDist)] .= 0
full_df.dS = normalize_yj(full_df.dS)

# calculate the VIFs by calculating the correlation matrix and taking the diagonal of its inverse
mat = Matrix(full_df[:,4:end])
cor_mat = cor(mat)
inv_cor_mat = inv(cor_mat)
vif_values = diag(inv_cor_mat)

# Create a heatmap of the correlation matrix, with values annotated
heatmap_trace = heatmap(
    z=cor_mat,
    x=names(full_df)[4:end],
    y=names(full_df)[4:end],
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
savefig(heatmap_fig, joinpath(data_dir, "correlation_heatmap_$(life_cycle).html"))

mlr_model = lm(
    @formula(
        dS ~
        AvgExpr +
        H3K27ac +
        H3K4me3 +
        H3K9me3 +
        ATAC +
        HasH3K27ac +
        HasH3K4me3 +
        HasH3K9me3 +
        HasATAC +
        SameChrom +
        ParalogDist
    ), 
    full_df)
mlr_table = DataFrame(coeftable(mlr_model))
rename!(mlr_table, [:Predictor, :Coef, :StdErr, :tvalue, :pvalue, :Lower95CI, :Upper95CI])
r_squared = r²(mlr_model)
valid_p_vals = findall(x -> !isnan(x), mlr_table.pvalue)
invalid_p_vals = setdiff(1:nrow(mlr_table), valid_p_vals)
mlr_table.pvalue_adj = fill(NaN, nrow(mlr_table))
mlr_table.pvalue_adj[valid_p_vals] = adjust(mlr_table[!, :pvalue][valid_p_vals], BenjaminiHochberg())

sort!(mlr_table, :pvalue_adj)
CSV.write(joinpath(data_dir, "mlr_results_$(life_cycle).csv"), mlr_table)
stats_file = open(joinpath(data_dir, "mlr_stats_$(life_cycle).txt"), "w")
write(stats_file, "r²: $r_squared\n\nVIF Values:\n")
for (i, col_name) in enumerate(names(full_df)[4:end])
    write(stats_file, "$col_name: $(vif_values[i])\n")
end
close(stats_file)