include("prelude.jl")

using .RepeatUtils

plot_save_dir = "../../dicty_data/saved_figures/"

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
addexpression!(ref_genome, CSV.read(expr_data_file, DataFrame))
paralog_df = CSV.read(paralog_file, DataFrame)
singleton_df = CSV.read(singleton_file, DataFrame)

# Convert TE genes to 'repeats' in the reference genome object:
te_genes_list = [gene for gene in ref_genome.genes[2] if gene.id in te_ids]
non_te_gene_list = [gene for gene in ref_genome.genes[2] if gene.id ∉ te_ids]
duplicates = [gene for gene in non_te_gene_list if gene.id in paralog_df.GeneID || gene.id in paralog_df.ParalogID]
singletons = [gene for gene in non_te_gene_list if gene.id in singleton_df.GeneID]

te_box = box(y = [isempty(gene.rnas) ? log2(0.5) : log2(mean(gene.rnas[1].expression) + 0.5) for gene in te_genes_list], name = "TE genes", color="blue")
dup_box = box(y = [isempty(gene.rnas) ? log2(0.5) : log2(mean(gene.rnas[1].expression) + 0.5) for gene in duplicates], name = "Paralogs")
sing_box = box(y = [isempty(gene.rnas) ? log2(0.5) : log2(mean(gene.rnas[1].expression) + 0.5) for gene in singletons], name = "Singletons")
non_te_box = box(y = [isempty(gene.rnas) ? log2(0.5) : log2(mean(gene.rnas[1].expression) + 0.5) for gene in non_te_gene_list], name = "Non-TE genes", color="red")

layout = Layout(yaxis=attr(title="log2(TPM + 0.5)", 
                           titlefont=attr(family = "Times New Roman", 
                           size=30),
                           gridcolor="lightgray",
                           gridwidth=1.5), 
                xaxis=attr(title="Gene type",
                           titlefont=attr(family = "Times New Roman", 
                           size=30),
                           gridcolor="lightgray",
                           gridwidth=1.5),
                plot_bgcolor="white",
                barmode="group",
                gridcolor="lightgray",
                gridwidth=1.5)
fig = plot([te_box, non_te_box], layout)

savefig(fig, joinpath(plot_save_dir, "gene_type_expression.html"))

non_te_expression = [mean(gene.rnas[1].expression) for gene in non_te_gene_list if length(gene.rnas) > 0]
te_expression = [mean(gene.rnas[1].expression) for gene in te_genes_list if length(gene.rnas) > 0]
MannWhitneyUTest(non_te_expression, te_expression)