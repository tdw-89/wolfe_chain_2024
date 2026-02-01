include("prelude.jl")

# Custom lib src:
using .ParalogUtils

# Paralog file
full_paralog_file = "../../dicty_data/AX4/genome_ver_2_7/biomart/genes_v52/d_disc_paralogs_biomart_ensembl_protist_ver_52.txt"

# Gene blacklist files
blacklist_file = "./blacklists/cds_blacklist_full.tsv"

# dictybase_cds_id_file = "../../dicty_data/AX4/genome_ver_2_7/fastas/dicty_primary_cds_ids.txt"
ensembl_cds_id_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/cds_ids.txt"

# Load paralog data
paralog_data = CSV.read(
    full_paralog_file, 
    DataFrame, 
    delim='\t'
    )
singleton_data = filter(row -> ismissing(row[15]), paralog_data)
singleton_ids = unique(singleton_data[!, 1])

# Save the unfiltered singleton list
CSV.write("../../dicty_data/singleton_list_unfilt.txt", DataFrame(GeneID = singleton_data[!,1]))

# Load blacklist
blacklist = CSV.read(blacklist_file, DataFrame)

# Load cds ID list
cds_ids = open(ensembl_cds_id_file) do file
    readlines(file)
end

###########
# Filter: #
###########

# Filter out gene splits, and genes with missing paralogy type annotation
filter!(row -> !ismissing(row[7]) && row[7] != "gene_split", paralog_data)
select!(paralog_data, [1,15,5,6,4,3,8])
rename!(paralog_data, ["GeneID", "ParalogID", "PercIDqt","PercIDtq","dS","dN", "LCA"])
filter!(row -> !ismissing(row.ParalogID), paralog_data)

# Filter out genes/paralogs in the blacklist or not in the CDS list
filter!(row -> row.GeneID in cds_ids && row.ParalogID in cds_ids, paralog_data)
filter!(row -> row.GeneID ∉ blacklist.GeneID && row.ParalogID ∉ blacklist.GeneID, paralog_data)
filter!(id -> id ∉ blacklist.GeneID, singleton_ids)

# Remove duplicates
unique!(paralog_data)

# Filter out genes/paralogs with a %ID below 30
filter!(row -> row[3] >= 30 && row[4] >= 30, paralog_data)

# Get RBH pairs
# Only use one round of 'rbh', otherwise you will include paralogs that have more recent duplicates,
# but are nonetheless paired with a more diverged duplicate (giving the false impression that their paired/averaged data
# represents the state of more diverged pairs that haven't experienced any duplication since their divergence)

rbh_df = rbh(paralog_data, scoring="mean")
rbh_id_pairs = collect(zip(rbh_df.GeneID, rbh_df.ParalogID))
filter!(row -> (row.GeneID, row.ParalogID) in rbh_id_pairs, paralog_data)

# NOTE: All rbh pairs are in the paralog data twice, once with the gene ID in the first column, and once with the paralog ID in the first column
# Proof:
#    rbh_id_paris_r = collect(zip(rbh_df.ParalogID, rbh_df.GeneID))
#    paralog_data_f = filter(row -> (row.GeneID, row.ParalogID) in rbh_id_pairs, paralog_data) #
#    paralog_data_r = filter(row -> (row.GeneID, row.ParalogID) in rbh_id_paris_r, paralog_data) #
#    sort!(paralog_data_r, :ParalogID)
#    sort!(paralog_data_f, :GeneID)
#    all(paralog_data_f.GeneID .== paralog_data_r.ParalogID)

# filter genes with no expression in any sample:
expr_data = "../../dicty_data/filtered/expr_data_filt_kallisto_ensembl52_single.tsv"
expr_data = CSV.read(expr_data, DataFrame)
filter!(row -> row.GeneID in expr_data.GeneID && row.ParalogID in expr_data.GeneID, paralog_data)
filter!(id -> id in expr_data.GeneID, singleton_ids)

CSV.write("../../dicty_data/filtered/paralog_filt.tsv", paralog_data, delim='\t', header=true)
CSV.write("../../dicty_data/filtered/singleton_filt.tsv", DataFrame(:GeneID => singleton_ids), delim='\t', header=true)