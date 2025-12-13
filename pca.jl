using MultivariateStats
using MultipleTesting
using CSV
using DataFrames
using CategoricalArrays
using PlotlyJS
using NaturalSort
using StatsBase

id_method = "dS" # "dS", "dN", or "id"
id_threshold = 3

# Custom lib src:
include("./custom_lib/load_gff.jl")
include("./custom_lib/genomic_data.jl")
include("./custom_lib/enrichment_utils.jl")
include("./custom_lib/misc_utils.jl")

# Genome data
gff_data = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"

# Expression data
expr_data_file = "../../dicty_data/filtered/expr_data_filt_kallisto_ensembl52_single.tsv"

# Peak files
chip_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_1_ensembl52/"
atac_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_2_ensembl52/"

# Paralog data
paralog_file = "../../dicty_data/filtered/paralog_filt.tsv"

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
expr_data.Avg = log.(expr_data.Avg .+ 0.5)

# Load peak data
peak_files = [filter(fn -> !contains(fn, r"_S[AB]_"), readdir(chip_peak_file_dir, join=true)); readdir(atac_peak_file_dir, join=true)]
filter!(fn -> endswith(fn, ".narrowPeak"), peak_files)
peak_data = binpeaks(peak_files, chrom_lengths_file)

# Load paralog
paralog_data = CSV.read(paralog_file, DataFrame)

# Filter the paralog data to just those that have expression data
filter!(row -> row.GeneID in expr_data.GeneID && row.ParalogID in expr_data.GeneID, paralog_data)

# Filter paralogs according to id threshold
if id_method in ["dS", "dN"]
    select!(paralog_data, ["GeneID", "ParalogID", id_method])
    filter!(row -> row[id_method] <= id_threshold, paralog_data)

elseif id_method == "id"
    mean_id = map(pair -> mean(pair), zip(paralog_data.PercIDqt, paralog_data.PercIDtq))
    insertcols!(paralog_data, 3, :MeanPerc => mean_id)
    filter!(row -> row.MeanPerc >= id_threshold, paralog_data)
    select!(paralog_data, ["GeneID", "ParalogID", "MeanPerc"])

else
    error("Invalid id_method: $id_method")

end

# Load reference genome object
ref_genome = loadgenome(gff_data, chrom_lengths_file)

# Add data to the reference
addexpression!(ref_genome, expr_data)
addtogenes!(ref_genome, peak_data)

# Create a dataframe for individual paralog data
full_df = DataFrame(GeneID = String.(collect(vcat(paralog_data.GeneID, paralog_data.ParalogID))), dS = vcat(paralog_data.dS, paralog_data.dS))

# Get the average expression for each gene
insertcols!(full_df, :AvgExpr => zeros(nrow(full_df)))
for (i, id) in enumerate(full_df.GeneID)
    full_df[i, :AvgExpr] = only(only(get(ref_genome, id).rnas).expression)

end

# Get the average enrichment for each gene for each mark

    # H3K27ac
insertcols!(full_df, :H3K27ac => zeros(nrow(full_df)))
for (i, id) in enumerate(full_df.GeneID)
    full_df[i, :H3K27ac] = mean([mean(getsiginrange(get(ref_genome, id), GeneRange(TSS(), TES(), 0, 0), ind)) for ind in [1, 4, 7]])

end

    # H3K4me3
insertcols!(full_df, :H3K4me3 => zeros(nrow(full_df)))
for (i, id) in enumerate(full_df.GeneID)
    full_df[i, :H3K4me3] = mean([mean(getsiginrange(get(ref_genome, id), GeneRange(TSS(), TES(), 0, 0), ind)) for ind in [2, 5, 8]])

end

    # H3K9me3
insertcols!(full_df, :H3K9me3 => zeros(nrow(full_df)))
for (i, id) in enumerate(full_df.GeneID)
    full_df[i, :H3K9me3] = mean([mean(getsiginrange(get(ref_genome, id), GeneRange(TSS(), TES(), 0, 0), ind)) for ind in [3, 6, 9]])

end

    # ATAC
insertcols!(full_df, :ATAC => zeros(nrow(full_df)))
for (i, id) in enumerate(full_df.GeneID)
    full_df[i, :ATAC] = mean([mean(getsiginrange(get(ref_genome, id), GeneRange(TSS(), TES(), 0, 0), ind)) for ind in [10, 11, 12]])

end


plot([box(y= full_df.AvgExpr[full_df.H3K9me3 .> 0], name="w/ H3K9me3"), box(y= full_df.AvgExpr[full_df.H3K9me3 .<= 0], name="w/o H3K9me3")], Layout(title="AvgExpr (w/ H3K9me3 vs w/o H3K9me3)"))
plot([box(y= full_df.AvgExpr[full_df.H3K27ac .> 0], name="w/ H3K27ac"), box(y= full_df.AvgExpr[full_df.H3K27ac .<= 0], name="w/o H3K27ac")], Layout(title="AvgExpr (w/ H3K27ac vs w/o H3K27ac)"))
plot([box(y= full_df.AvgExpr[full_df.H3K4me3 .> 0], name="w/ H3K4me3"), box(y= full_df.AvgExpr[full_df.H3K4me3 .<= 0], name="w/o H3K4me3")], Layout(title="AvgExpr (w/ H3K4me3 vs w/o H3K4me3)"))
plot([box(y= full_df.AvgExpr[full_df.ATAC .> 0], name="w/ ATAC"), box(y= full_df.AvgExpr[full_df.ATAC .<= 0], name="w/o ATAC")], Layout(title="AvgExpr (w/ ATAC vs w/o ATAC)"))

qs = levelcode.(cut(full_df.dS, 10))
plot([box(y= full_df[qs .<= 2,:].AvgExpr[full_df[qs .<= 2,:].H3K9me3 .> 0], name="w/ H3K9me3"), box(y= full_df[qs .<= 2,:].AvgExpr[full_df[qs .<= 2,:].H3K9me3 .<= 0], name="w/o H3K9me3")], Layout(title="AvgExpr (w/ H3K9me3 vs w/o H3K9me3, low-dS)"))

# PCA ANALYSIS: #

# Z-score normalize the data (except dS)
full_df.AvgExpr = zscore(full_df.AvgExpr)
full_df.H3K27ac = zscore(full_df.H3K27ac)
full_df.H3K4me3 = zscore(full_df.H3K4me3)
full_df.H3K9me3 = zscore(full_df.H3K9me3)
full_df.ATAC = zscore(full_df.ATAC)

# # PCA

input_data = Matrix(full_df[qs .<= 2,:][:,3:end])'
pca_model = fit(PCA, input_data, maxoutdim=2)
pca_data = MultivariateStats.predict(pca_model, input_data)

# # Plot
fig = plot(
    scatter(x=pca_data[1, :], y=pca_data[2, :], mode="markers", text=full_df.GeneID[qs .<= 2], marker_color=full_df.AvgExpr[qs .<= 2], marker_colorscale="Viridis", marker_colorbar_title="AvgExpr"),
    Layout(title="PCA of expression and enrichment data (low dS)", xaxis_title="PC1", yaxis_title="PC2")
)