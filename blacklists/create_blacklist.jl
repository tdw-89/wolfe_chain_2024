using CSV
using DataFrames

cd("blacklists")

pseudogene_list = "ensembl_pseudogenes.tsv"
te_list_file = "ensembl_te_ids_with_predicted.tsv"
ax4_dup_file = "ax4_duplication_genes.tsv"

# gene_info_df = CSV.read(db_gene_info, DataFrame)
# rename!(gene_info_df, ["GeneID", "GeneName", "Synonyms", "GeneProducts"])
# db_pseudo_genes = filter(row -> contains(row.GeneName, "_ps") || (!ismissing(row.GeneProducts) && contains(lowercase(row.GeneProducts), "pseudogene")), gene_info_df).GeneID
pseudo_genes_df = CSV.read(pseudogene_list, DataFrame)
te_list = CSV.read(te_list_file, DataFrame)
# ax4_dup = CSV.read(ax4_dup_file, DataFrame, delim=',')

CSV.write("cds_blacklist_full.tsv", DataFrame(GeneID=unique(vcat(pseudo_genes_df.GeneID, te_list.GeneID, ax4_dup.GeneID))))