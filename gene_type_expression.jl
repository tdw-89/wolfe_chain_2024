include("prelude.jl")

using .RepeatUtils

# Paralog data
paralog_file = "../../dicty_data/filtered/paralog_filt.tsv"
singleton_file = "../../dicty_data/filtered/singleton_filt.tsv"

# TE files
te_id_file = "blacklists/ensembl_te_ids_with_predicted.tsv"
te_ids = open(te_id_file) do file
    readlines(file)
end

# Expression data
expr_data_file = "../../dicty_data/filtered/expr_data_filt_kallisto_ensembl52_single_with_TEs.tsv"

# Genome data
gff_source = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"

# Load genome data
ref_genome = loadgenome(gff_source, chrom_lengths_file)
expr_data = CSV.read(expr_data_file, DataFrame)
rename!(expr_data, ["GeneID", "V", "S", "M", "F"])
select!(expr_data, ["GeneID", "V", "M", "F"])  # stages with ChIP-seq data
addexpression!(ref_genome, expr_data)
paralog_df = CSV.read(paralog_file, DataFrame)
singleton_df = CSV.read(singleton_file, DataFrame)

# Convert TE genes to 'repeats' in the reference genome object:
te_genes_list = [gene for gene in ref_genome.genes[2] if gene.id in te_ids]
non_te_gene_list = [gene for gene in ref_genome.genes[2] if gene.id ∉ te_ids]
duplicates = [gene for gene in non_te_gene_list if gene.id in paralog_df.GeneID || gene.id in paralog_df.ParalogID]
singletons = [gene for gene in non_te_gene_list if gene.id in singleton_df.GeneID]

non_te_expression = [mean(gene.rnas[1].expression) for gene in non_te_gene_list if length(gene.rnas) > 0]
te_expression = [mean(gene.rnas[1].expression) for gene in te_genes_list if length(gene.rnas) > 0]

# Supplementary figure S2: expression of non-TE vs TE genes, log2(mean TPM over V/M/F + 0.5)
write_stats("te_gene_expression", DataFrame([mwu_row("Non-TE genes vs TE genes, log2(TPM + 0.5)",
    log2.(non_te_expression .+ 0.5), log2.(te_expression .+ 0.5))]))

fig = plot([fig_box(log2.(te_expression .+ 0.5), "TE genes", BLUE), fig_box(log2.(non_te_expression .+ 0.5), "Non-TE genes", RED)],
           fig_layout(yaxis=fig_axis(title="log<sub>2</sub>(TPM + 0.5)"), xaxis=fig_axis(showgrid=false)))
save_figure(fig, "figS02"; width=700, height=500)
