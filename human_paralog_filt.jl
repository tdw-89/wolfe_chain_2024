using CSV
using DataFrames

include("custom_lib/load_gff.jl")
include("custom_lib/paralog_utils.jl")

# Full paralog info file
human_paralog_info = "../../../../data/mammals/primates/h_sapiens/Ensembl_99/paralog_info_GRCh38_ensembl_99.txt"
cds_id_file = "../../../../data/mammals/primates/h_sapiens/Ensembl_99/CDS_IDs.txt"
nmd_id_file = "../../../../data/mammals/primates/h_sapiens/Ensembl_99/nmd_candidate_ids.txt"

# Load the paralog data
paralog_data = CSV.read(human_paralog_info, DataFrame)
cds_df = CSV.read(cds_id_file, DataFrame, header=false)
nmd_df = CSV.read(nmd_id_file, DataFrame, header=false)
select!(paralog_data, [1,11,7,8,5,6])
unique!(paralog_data)
rename!(paralog_data, ["GeneID", "ParalogID", "PercIDqt", "PercIDtq", "dS","dN"])

# Identify the singletons
singletons = unique(select(filter(row -> ismissing(row.ParalogID), paralog_data), :GeneID))
filter!(row -> row.GeneID ∈ cds_df[!,1] && row.GeneID ∉ nmd_df[!,1], singletons)
CSV.write("data/filtered/human_singletons_filt.csv", singletons)

# Filter the paralog data
filter!(row -> !ismissing(row[3]) && !ismissing(row[4]) && !ismissing(row.dS), paralog_data) # Filter to those with %ID info
filter!(row -> min(row[3],row[4]) >= 30, paralog_data) # Filter to those with %ID >= 30
filter!(row -> row.GeneID ∈ cds_df[!,1] && row.ParalogID ∈ cds_df[!,1], paralog_data) # Filter to those that are in the coding sequence file
filter!(row -> row.GeneID ∉ nmd_df[!,1] && row.ParalogID ∉ nmd_df[!,1], paralog_data) # Filter to those that are not in the NMD candidate file

# Filter to the reciprocal best hits
rbh_df = rbh(paralog_data, scoring="mean")
rbh_id_pairs = collect(zip(rbh_df.GeneID, rbh_df.ParalogID))
filter!(row -> (row.GeneID, row.ParalogID) in rbh_id_pairs, paralog_data)

CSV.write("data/filtered/human_paralog_info_filt.csv", paralog_data)