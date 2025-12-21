include("prelude.jl")

# Custom lib src:
include("./custom_lib/paralog_utils.jl")

# Filtering parameters
filter_not_expressed = true
ds_est = false
for_ds_calc = false
yn00 = false
ng86 = false
lwl85 = false

# Paralog file
full_paralog_file = "../../dicty_data/AX4/genome_ver_2_7/biomart/genes_v52/d_disc_paralogs_biomart_ensembl_protist_ver_52.txt"

# Gene blacklist files
blacklist_file = "./blacklists/cds_blacklist_full.tsv"
# dictybase_cds_id_file = "../../dicty_data/AX4/genome_ver_2_7/fastas/dicty_primary_cds_ids.txt"
ensembl_cds_id_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/cds_ids.txt"

# Load paralog data
paralog_data = CSV.read(full_paralog_file, DataFrame, delim='\t')
singleton_data = filter(row -> ismissing(row[15]), paralog_data)
singleton_ids = unique(singleton_data[!, 1])
ds_estimates = CSV.read("../../dicty_data/filtered/dS_df.csv", DataFrame)

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

# # OLD METHOD:
# input_df = copy(paralog_data)
# rbh_df = rbh(input_df, scoring="mean")
# rbh_ids = vcat(rbh_df.GeneID, rbh_df.ParalogID)
# filter!(row -> row.GeneID ∉ rbh_ids && row.ParalogID ∉ rbh_ids, input_df)

# while nrow(input_df) > 0 #?
#     append!(rbh_df, rbh(input_df, scoring="mean"))
#     rbh_ids = vcat(rbh_df.GeneID, rbh_df.ParalogID)
#     filter!(row -> row.GeneID ∉ rbh_ids && row.ParalogID ∉ rbh_ids, input_df)
    
# end
# rbh_id_pairs = collect(zip(rbh_df.GeneID, rbh_df.ParalogID))
# filter!(row -> (row.GeneID, row.ParalogID) in rbh_id_pairs, paralog_data) #

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


# Write filtered paralog data
if ds_est

    if for_ds_calc
        # Write the pair IDs for the dS calculation pipeline
        id_pairs = select(paralog_data, [1,2])
        CSV.write("../../dicty_data/filtered/paralog_filt_id_pairs.tsv", id_pairs, delim='\t', header=false)

    else

        method_ind = [3,5,7][only(findall([yn00, ng86, lwl85]))]
        ds_estimate_sets = [Set([ds_estimates[i, 1], ds_estimates[i, 2]]) for i in 1:nrow(ds_estimates)]

        for i in 1:nrow(paralog_data)
            temp_set = Set([paralog_data[i, 1], paralog_data[i, 2]])
            ds_ind = findfirst(set -> set == temp_set, ds_estimate_sets)
            paralog_data.dS[i] = ds_estimates[ds_ind, method_ind]

        end

        CSV.write("../../dicty_data/filtered/paralog_filt.tsv", paralog_data, delim='\t', header=true)
    end
else
    CSV.write("../../dicty_data/filtered/paralog_filt.tsv", paralog_data, delim='\t', header=true)

end

CSV.write("../../dicty_data/filtered/singleton_filt.tsv", DataFrame(:GeneID => singleton_ids), delim='\t', header=true)

#=
OLD_METHOD:

    # original filtering parameters:
    ds_threshold = 2.5
    filter_tes = true
    filter_ncRNAs = true
    filter_all_nc = true
    filter_plasmids = true
    filter_pseudogenes = true
    filter_ax4_dup_region = true
    filter_mito = true
    filter_not_expressed = false

    select!(paralog_data, [1,15,5,6,4])
    rename!(paralog_data, ["GeneID", "ParalogID", "PercIDqt","PercIDtq","dS"])
    filter!(row -> !ismissing(row.ParalogID), paralog_data)
    unique!(paralog_data)
    filter!(row -> !ismissing(row.dS), paralog_data)
    filter!(row -> row.dS <= ds_threshold, paralog_data)

    rbh_df = rbh(paralog_data, scoring="mean")
    filter!(row -> row.mean_perc >= 30, rbh_df)
    rbh_id_pairs = collect(zip(rbh_df.GeneID, rbh_df.ParalogID))
    filter!(row -> (row.GeneID, row.ParalogID) in rbh_id_pairs, paralog_data)

NEWOLD_METHOD:

    paralog_data = CSV.read("../../dicty_data/wang_et_al/paralogs/Dicty.GeneParalogPairs.v51.csv", DataFrame)
    filter!(row -> row.avg_pct_id >= 30, paralog_data)
    filter!(row -> row.ddiscoideum_eg_paralog_orthology_type != "gene_split", paralog_data)
    select!(paralog_data, [2,4,10,11,13])
    rename!(paralog_data, ["GeneID", "ParalogID", "PercIDqt","PercIDtq","dS"])
    paralog_data.dS = parse.(Float64, paralog_data.dS)
    unique!(paralog_data)
    filter!(row -> row.dS <= 5, paralog_data)
=#