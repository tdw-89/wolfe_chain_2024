using Pkg
Pkg.activate("BioinfoTools/")
using BioinfoTools.LoadGFF
using BioinfoTools.EnrichmentUtils
using BioinfoTools.GenomicData
using BioinfoTools.GenomeTypes
import BioinfoTools.MiscUtils: normalize_yj

using CategoricalArrays
using CSV
using DataFrames
using GLM
using MultipleTesting
using PlotlyJS
using StatsBase
using Serialization

data_dir = "../../dicty_data/"
reload_peak_data = false

chip_peak_file_dir = joinpath(data_dir, "mammals/primates/h_sapiens/ENCODE_histone_mods/")
gff_source = joinpath(data_dir, "mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.99.gff3")
chrom_lengths_file = joinpath(data_dir, "mammals/primates/h_sapiens/Ensembl_99/chromosome_lengths.txt")
paralog_file = joinpath(data_dir, "filtered/human_paralog_info_filt.csv")
cds_id_file = joinpath(data_dir, "mammals/primates/h_sapiens/Ensembl_99/CDS_IDs.txt")
nmd_id_file = joinpath(data_dir, "mammals/primates/h_sapiens/Ensembl_99/nmd_candidate_ids.txt")

upstream, downstream = -500, 500
id_threshold = 3
zlog_normalize(u::Vector{Float64}) = zscore(log.(u))

# Load the genome data:
ref_genome = loadgenome(gff_source, chrom_lengths_file)

# Load the paralog data:
paralog_data = CSV.read(paralog_file, DataFrame)
cds_df = CSV.read(cds_id_file, DataFrame, header=false)
nmd_df = CSV.read(nmd_id_file, DataFrame, header=false)

# Filter the paralogs:
select!(paralog_data, ["GeneID", "ParalogID", "dS"])
filter!(row ->
    row.dS <= id_threshold &&
    row.GeneID ∈ cds_df.Column1 && 
    row.ParalogID ∉ nmd_df.Column1, 
    paralog_data
    )

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
idxs = collect(1:length(ref_genome.genes[2][1].samples))
h3k9me3_means = Vector{Float64}()

for i in 1:nrow(paralog_data)
    gid = paralog_data.GeneID[i]
    pid = paralog_data.ParalogID[i]
    gid_sig_h3k9me3 = mean([mean(getsiginrange(get(ref_genome, gid), GeneRange(TSS(), TES(), upstream, downstream), ind)) for ind in idxs])
    pid_sig_h3k9me3 = mean([mean(getsiginrange(get(ref_genome, pid), GeneRange(TSS(), TES(), upstream, downstream), ind)) for ind in idxs])
    mean_sig = mean([gid_sig_h3k9me3, pid_sig_h3k9me3])
    push!(h3k9me3_means, mean_sig)
end

# normalize the data
paralog_data.H3K9me3 = normalize_yj(h3k9me3_means)
paralog_data.dS = normalize_yj(paralog_data.dS)

# Fit a linear model
mlr_formula = @formula(dS ~ H3K9me3)
mlr_model = lm(mlr_formula, paralog_data)
mlr_table = DataFrame(coeftable(mlr_model))
rename!(mlr_table, [:Predictor, :Coef, :StdErr, :tvalue, :pvalue, :Lower95CI, :Upper95CI])
r_squared = r²(mlr_model)

CSV.write(joinpath(data_dir, "human_paralog_h3k9me3_mlr_results.csv"), mlr_table)
open(joinpath(data_dir, "human_paralog_h3k9me3_mlr_r2.txt"), "w") do io
    write(io, "R²: $(r_squared)\n")
end