include("prelude.jl")

using CategoricalArrays

LOW_DS_QUANT_NUM_CUTOFF = 2
global_range = GeneRange(TSS(), TES(), -500, 500)

using .EnrichmentUtils
using .MiscUtils

function in_significant_range(gene::Gene, sig_range, ind)
    return sum(getsiginrange(gene, sig_range, ind)) > 0
end

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
expr_data[!,2:end] = log2.(expr_data[!,2:end] .+ 0.5)

# Load peak data
peak_files = [filter(fn -> !contains(fn, r"_S[AB]_"), readdir(chip_peak_file_dir, join=true)); readdir(atac_peak_file_dir, join=true)]
filter!(fn -> endswith(fn, ".narrowPeak"), peak_files)
peak_data = binpeaks(peak_files, chrom_lengths_file)

# Load paralog data
paralog_data = CSV.read(paralog_file, DataFrame)

# Filter paralogs according to id threshold
select!(paralog_data, ["GeneID", "ParalogID", "dS"])
filter!(row -> row.dS <= 3.0, paralog_data)

# Filter to pairs which have expression data for both
filter!(row -> row.GeneID in expr_data.GeneID && row.ParalogID in expr_data.GeneID, paralog_data) 

# Load reference genome object
ref_genome = loadgenome(gff_data, chrom_lengths_file)

# Add data to reference
addexpression!(ref_genome, expr_data)
addtogenes!(ref_genome, peak_data)

# Get the dS values and quantiles for the filtered pairs
quantile_labels = cut(paralog_data[!, 3], 10)
quantile_vals = levelcode.(quantile_labels)

# Individual expression vs K9me3 enrichment
    # Low dS paralog genes
paralog_genes_low = 
    get(ref_genome, String.(vcat(
        paralog_data.GeneID[quantile_vals .<= LOW_DS_QUANT_NUM_CUTOFF], 
        paralog_data.ParalogID[quantile_vals .<= LOW_DS_QUANT_NUM_CUTOFF]
        )))
paralog_gene_expr_low = [mean(gene.rnas[1].expression) for gene in paralog_genes_low]
paralog_gene_has_k9me3_low_1 = [in_significant_range(gene, global_range, 3) for gene in paralog_genes_low]
paralog_gene_has_k9me3_low_2 = [in_significant_range(gene, global_range, 6) for gene in paralog_genes_low]
paralog_gene_has_k9me3_low_3 = [in_significant_range(gene, global_range, 9) for gene in paralog_genes_low]
paralog_gene_has_k9me3_low = paralog_gene_has_k9me3_low_1 .| paralog_gene_has_k9me3_low_2 .| paralog_gene_has_k9me3_low_3

# Table for plotting (in R):

low_ds_df = DataFrame(
    :GeneID => [gene.id for gene in paralog_genes_low],
    :HasK9me3 => paralog_gene_has_k9me3_low,
    :HasK27ac => false,
    :HasK4me3 => false,
    :HasATAC => false,
    :Expr => paralog_gene_expr_low
)

# Add H3K4me3, H3K27ac, and ATAC data
k27ac_inds = [1,4,7]
k4me3_inds = [2,5,8]
atac_inds = [10,11,12]

Threads.@threads for i in eachindex(paralog_genes_low)

    gene = paralog_genes_low[i]
    low_ds_df.HasK27ac[i] = sum([sum(getsiginrange(gene, global_range, ind)) for ind in k27ac_inds]) > 0
    low_ds_df.HasK4me3[i] = sum([sum(getsiginrange(gene, global_range, ind)) for ind in k4me3_inds]) > 0
    low_ds_df.HasATAC[i] = sum([sum(getsiginrange(gene, global_range, ind)) for ind in atac_inds]) > 0
end 

CSV.write("../../dicty_data/low_ds_df.csv", low_ds_df)