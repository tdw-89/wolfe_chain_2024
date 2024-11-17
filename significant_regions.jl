using GLM
using CSV
using DataFrames
using HypothesisTests
using MultipleTesting
using PlotlyJS

# Custom lib src:
include("./custom_lib/load_gff.jl")
include("./custom_lib/genomic_data.jl")
include("./custom_lib/enrichment_utils.jl")

const flanking_region = 500

# Functions:
function full_region(gene, sample_ind, flank_length)
    signal = [getsiginrange(gene, GeneRange(TSS(), TSS(), -flank_length, -1), sample_ind); 
              to_percent(Float64.(getsiginrange(gene, GeneRange(TSS(), TES()), sample_ind))); 
              getsiginrange(gene, GeneRange(TES(), TES(), 1, flank_length), sample_ind)]
    return signal
end

function scan_tss(df::DataFrame)

    tss_ind = findfirst(df.Pos .== "1%")
    
    if isnothing(tss_ind)
        error("TSS position ('1%') not found in dataframe")
    
    end

    # Scan upstream:
    upstream_inds = []
    upstream_pos = []
    for i in tss_ind:-1:1
        
        if df.Pval[i] <= 0.05
            push!(upstream_inds, i)
            push!(upstream_pos, df.Pos[i])

        else
            break
        end
    end

    # Scan downstream:
    downstream_inds = []
    downstream_pos = []
    for i in tss_ind+1:length(df.Pos)

        if df.Pval[i] <= 0.05
            push!(downstream_inds, i)
            push!(downstream_pos, df.Pos[i])

        else
            break
        end
    end

    return df[sort([upstream_inds; downstream_inds]),:]
end

function sig_regions(gene_list, expr_vec, sample_inds, plot_label="Unknown"; show_plot=true)

    sig_mats_body = []
    sig_mats_upstream = []
    sig_mats_downstream = []
    for ind in sample_inds
        sig_mat_body = zeros(length(gene_list), 100)
        sig_mat_upstream = zeros(length(gene_list), flanking_region)
        sig_mat_downstream = zeros(length(gene_list), flanking_region)

        for i in axes(sig_mat_body, 1)
            sig = to_percent(Float64.(getsiginrange(gene_list[i], GeneRange(TSS(), TES(), 0, 0), ind)))
            sig_mat_body[i, :] = sig

        end

        for i in axes(sig_mat_upstream, 1)
            sig = Float64.(getsiginrange(gene_list[i], GeneRange(TSS(), TSS(), -flanking_region, -1), ind))
            sig_mat_upstream[i, :] = sig

        end

        for i in axes(sig_mat_downstream, 1)
            sig = Float64.(getsiginrange(gene_list[i], GeneRange(TES(), TES(), 1, flanking_region), ind))
            sig_mat_downstream[i, :] = sig

        end

        push!(sig_mats_body, sig_mat_body)
        push!(sig_mats_upstream, sig_mat_upstream)
        push!(sig_mats_downstream, sig_mat_downstream)

    end

    # Average the matrices
    sig_mat_body = mean(sig_mats_body)
    sig_mat_upstream = mean(sig_mats_upstream)
    sig_mat_downstream = mean(sig_mats_downstream)

    # return sig_mat_body, sig_mat_upstream, sig_mat_downstream

    body_df = DataFrame(Pos = string.(collect(1:100)) .* "%", Cor = zeros(100), Pval = zeros(100))
    for j in axes(sig_mat_body, 2)
        lin_cor = cor(expr_vec, sig_mat_body[:, j])
        cor_pval = pvalue(CorrelationTest(expr_vec, sig_mat_body[:, j]))
        body_df[j, :Cor] = lin_cor
        body_df[j, :Pval] = cor_pval

    end

    upstream_df = DataFrame(Pos = string.(collect(-flanking_region:-1)), Cor = zeros(flanking_region), Pval = zeros(flanking_region))
    for j in axes(sig_mat_upstream, 2)
        lin_cor = cor(expr_vec, sig_mat_upstream[:, j])
        cor_pval = pvalue(CorrelationTest(expr_vec, sig_mat_upstream[:, j]))
        upstream_df[j, :Cor] = lin_cor
        upstream_df[j, :Pval] = cor_pval

    end

    downstream_df = DataFrame(Pos = "+" .* string.(collect(1:flanking_region)), Cor = zeros(flanking_region), Pval = zeros(flanking_region))
    for j in axes(sig_mat_downstream, 2)
        lin_cor = cor(expr_vec, sig_mat_downstream[:, j])
        cor_pval = pvalue(CorrelationTest(expr_vec, sig_mat_downstream[:, j]))
        downstream_df[j, :Cor] = lin_cor
        downstream_df[j, :Pval] = cor_pval

    end

    full_df = vcat(upstream_df, body_df, downstream_df)
    full_df.Pval = adjust(full_df.Pval, BenjaminiHochberg())
    # full_df.Pval = map(pval -> min(10, -log10(pval)), full_df.Pval)

    if show_plot
        display(plot([scatter(x=1:nrow(full_df), y=full_df.Cor, name="Correlation"), 
                    scatter(x=1:nrow(full_df), y=full_df.Pval, name="p-value (adj.)"),
                    # scatter(x=1:nrow(full_df), y=fill(1.3010299956639813, nrow(full_df)), name="Significance cutoff")], Layout(
                    scatter(x=1:nrow(full_df), y=fill(0.05, nrow(full_df)), name="Significance cutoff")], 
                    Layout(
            title=plot_label,
            xaxis=attr(tickmode="array", 
                       ticktext=vcat(upstream_df.Pos, body_df.Pos, downstream_df.Pos)[collect(100:100:nrow(full_df))], 
                       tickvals=collect(100:100:nrow(full_df)),
                       gridcolor="lightgray",
                       gridwidth=1.5),
            yaxis=attr(gridcolor="lightgray",
                       gridwidth=1.5),
            plot_bgcolor="white"
            )))

    end
    return scan_tss(full_df)
end

# Peak files
chip_peak_file_dir = "../../../../data/wang_et_al/processed/run_1_ensembl52/"
# chip_peak_file_dir = "../../../../data/wang_et_al/processed/run_16/"
# h3_peak_dir = "../../../../data/wang_et_al/processed/run_9/peak_files_H3"
atac_peak_file_dir = "../../../../data/wang_et_al/processed/run_2_ensembl52/"
# peak_dir = "../../../../data/zebrafish/data_from_freddy/"

# Genome data
# ensembl 52 genome data
gff_data = "../../../../data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../../../data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
ensembl_cds_id_file = "../../../../data/AX4/genome_ver_2_7/ensembl_52/cds_ids.txt"
blacklist_file = "./data/blacklist_with_tes.csv"
te_dist_file = "./data/te_distance.csv"

# zebrafish genome data
# gff_file = "../../../../data/zebrafish/Danio_rerio.GRCz10.91.gff3"
# chrom_lengths_file = "../../../../data/zebrafish/chromosome_lengths_z10.txt"

# Expression data
expr_data_file = "./data/filtered/expr_data_filt_kallisto_ensembl52_single.tsv"

# MAIN SCRIPT: #

# Load genome data
ref_genome = loadgenome(gff_data, chrom_lengths_file)
# ref_genome = loadgenome(gff_file, chrom_lengths_file)

# Load cds ID list
cds_ids = open(ensembl_cds_id_file) do file
    readlines(file)
end

# Load blacklist
blacklist = CSV.read(blacklist_file, DataFrame)

# Load peak data
peak_files = [filter(fn -> !contains(fn, r"_S[AB]_"), readdir(chip_peak_file_dir, join=true)); readdir(atac_peak_file_dir, join=true)]
filter!(fn -> endswith(fn, ".narrowPeak"), peak_files)
peak_data = binpeaks(peak_files, chrom_lengths_file)
addtogenes!(ref_genome, peak_data)

# Load expression data
expr_data = CSV.read(expr_data_file, DataFrame)

# Rename the expression columns to a letter indicating with life-cycle stage the sample
# was taken from.
rename!(expr_data, ["GeneID", "V", "S", "M", "F"])

# Remove 'S' stage expression data, and rearrange so that the samples are in the order they
# will be for the peak data
# select!(expr_data, ["GeneID", "F", "M", "V"])

# Average gene expression counts across all life-cycle stages
insertcols!(expr_data, :Avg => mean.(eachrow(expr_data[:, 2:end])))
select!(expr_data, ["GeneID", "Avg"])

# log-transform data with an added pseudocount of 0.5
expr_data.Avg = log.(expr_data.Avg .+ 0.5)

# Add expression data to genome object
addexpression!(ref_genome, expr_data)

# Filter to genes that have expression data and have 500 bp upstream and 500 downstream before a chromosome end
filtered_gene_list = [gene for gene in ref_genome.genes[2] if has_expr(gene) && siginrange(gene, GeneRange(TSS(), TES(), -flanking_region, flanking_region), peak_data=true)]

expr_vec = [only(gene.rnas[1].expression) for gene in filtered_gene_list]
k27_inds = [1,4,7]
k4_inds = [2,5,8]
k9me3_inds = [3,6,9]
atac_inds = [10,11,12]

# Check the distribution among the top 10% and bottom 10% of genes
expr_vals = [only(only(gene.rnas).expression) for gene in filtered_gene_list]
qs = levelcode.(cut(expr_vals, 10))
top_ten_percent = filtered_gene_list[findall(qs .== 10)]
bottom_ten_percent = filtered_gene_list[findall(qs .== 1)]
sigs_top_10_k4 = [mean([full_region(gene, ind, flanking_region) for ind in k4_inds]) for gene in top_ten_percent]
sigs_bottom_10_k4 = [mean([full_region(gene, ind, flanking_region) for ind in k4_inds]) for gene in bottom_ten_percent]
display(plot([scatter(x=1:length(sigs_top_10_k4), y=mean(sigs_top_10_k4), name="Top 10%"), scatter(x=1:length(sigs_bottom_10_k4), y=mean(sigs_bottom_10_k4), name="Bottom 10%")], Layout(title="K4me3")))
sigs_top_10_k27 = [mean([full_region(gene, ind, flanking_region) for ind in k27_inds]) for gene in top_ten_percent]
sigs_bottom_10_k27 = [mean([full_region(gene, ind, flanking_region) for ind in k27_inds]) for gene in bottom_ten_percent]
display(plot([scatter(x=1:length(sigs_top_10_k27), y=mean(sigs_top_10_k27), name="Top 10%"), scatter(x=1:length(sigs_bottom_10_k27), y=mean(sigs_bottom_10_k27), name="Bottom 10%")], Layout(title="K27ac")))
sigs_top_10_k9 = [mean([full_region(gene, ind, flanking_region) for ind in k9me3_inds]) for gene in top_ten_percent]
sigs_bottom_10_k9 = [mean([full_region(gene, ind, flanking_region) for ind in k9me3_inds]) for gene in bottom_ten_percent]
display(plot([scatter(x=1:length(sigs_top_10_k9), y=mean(sigs_top_10_k9), name="Top 10%"), scatter(x=1:length(sigs_bottom_10_k9), y=mean(sigs_bottom_10_k9), name="Bottom 10%")], Layout(title="K9me3")))
sigs_top_10_atac = [mean([full_region(gene, ind, flanking_region) for ind in atac_inds]) for gene in top_ten_percent]
sigs_bottom_10_atac = [mean([full_region(gene, ind, flanking_region) for ind in atac_inds]) for gene in bottom_ten_percent]
display(plot([scatter(x=1:length(sigs_top_10_atac), y=mean(sigs_top_10_atac), name="Top 10%"), scatter(x=1:length(sigs_bottom_10_atac), y=mean(sigs_bottom_10_atac), name="Bottom 10%")], Layout(title="ATAC")))

te_dist_df = CSV.read(te_dist_file, DataFrame)

k27_df = sig_regions(filtered_gene_list, expr_vec, k27_inds, "K27ac")
k4_df = sig_regions(filtered_gene_list, expr_vec, k4_inds, "K4me3")
k9_df = sig_regions(filtered_gene_list, expr_vec, k9me3_inds, "K9me3")
atac_df = sig_regions(filtered_gene_list, expr_vec, atac_inds, "ATAC")

sig_region_df = DataFrame("Mark" => ["K27ac", "K4me3", "K9me3", "ATAC"],
                        "Start" => [k27_df.Pos[1], k4_df.Pos[1], k9_df.Pos[1], atac_df.Pos[1]],
                        "End" => [k27_df.Pos[end], k4_df.Pos[end], k9_df.Pos[end], atac_df.Pos[end]])

CSV.write("./data/sig_regions.csv", sig_region_df)

# plot([scatter(x=1:length(k27_df.Cor),y=k27_df.Cor), scatter(x=1:length(k27_df.Pval), y=k27_df.Pval)], Layout(title="K27ac"))
# plot([scatter(x=1:length(k4_df.Cor),y=k4_df.Cor), scatter(x=1:length(k4_df.Pval), y=k4_df.Pval)], Layout(title="K4me3"))
# plot([scatter(x=1:length(k9_df.Cor),y=k9_df.Cor), scatter(x=1:length(k9_df.Pval), y=k9_df.Pval)], Layout(title="K9me3"))
# plot([scatter(x=1:length(atac_df.Cor),y=atac_df.Cor), scatter(x=1:length(atac_df.Pval), y=atac_df.Pval)], Layout(title="ATAC"))