# Run first
using CSV
using DataFrames

# Filtering parameters
filter_zeros = true
keep_tes = true

# Gene blacklist files
blacklist_file = "./blacklists/cds_blacklist_full.tsv"
ensembl_cds_id_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/cds_ids.txt"

# Expression data
# expr_data_file = "../../dicty_data/wang_et_al/RNA/kallisto_output_2/normalized_tpm_gene_names_old.tsv"
expr_data_file = "../../dicty_data/wang_et_al/RNA/kallisto_output_ensembl52/normalized_tpm_gene_names.tsv"

# Load expression data
expr_data = CSV.read(expr_data_file, DataFrame)

# Load blacklist
blacklist = CSV.read(blacklist_file, DataFrame)
blacklist = blacklist.GeneID

# Load cds ID list
cds_ids = open(ensembl_cds_id_file) do file
    readlines(file)
end

# If keeping TEs, load the TE ID list
if keep_tes
    te_id_file = "blacklists/ensembl_te_ids_with_predicted.tsv"
    te_ids = open(te_id_file) do file
        readlines(file)
    end

    blacklist = filter(id -> id ∉ te_ids, blacklist)
end

# Combine transcript counts from the same gene
id_col = "GeneID" in names(expr_data) ? "GeneID" : "gene"
filter!(row -> !ismissing(row[id_col]), expr_data)
col_types = map(col -> typeof(col[1]), eachcol(expr_data))
first_data_col = findfirst(col_types .== Float64)
unique_genes = unique(expr_data[!, id_col])
expr_data_combined = DataFrame(hcat(unique_genes, zeros(length(unique_genes), ncol(expr_data) - (first_data_col - 1))), :auto)
rename!(expr_data_combined, ["GeneID", names(expr_data[!, first_data_col:end])...])
for (i, gene) in enumerate(unique_genes)
    gene_inds = findall(expr_data[!, id_col] .== gene)
    expr_data_combined[i, 2:end] = sum(Matrix(expr_data[gene_inds, first_data_col:end]), dims=1)

end

# Filter

filter!(row -> row.GeneID ∉ blacklist, expr_data_combined)
filter!(row -> row.GeneID in cds_ids, expr_data_combined)

if filter_zeros
    filter!(row -> sum(row[2:end]) > 0, expr_data_combined)

end

# Average replicates:
rep_pairs = [(i, i + 1) for i in 2:2:ncol(expr_data_combined)]
for (i, pair) in enumerate(rep_pairs)
    expr_data_combined[!, pair[1]] = (expr_data_combined[!, pair[1]] .+ expr_data_combined[!, pair[2]]) ./ 2

end

select!(expr_data_combined, Not([pair[2] for pair in rep_pairs]))

CSV.write("../../dicty_data/filtered/expr_data_filt_kallisto_ensembl52_single$(keep_tes ? "_with_TEs" : "").tsv", expr_data_combined, delim='\t')