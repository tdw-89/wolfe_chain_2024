using CSV
using DataFrames

# Filtering parameters
filter_tes = true
filter_ncRNAs = true
filter_all_nc = true
filter_plasmids = true
filter_pseudogenes = true
filter_ax4_dup_region = true
filter_mito = true

# Blacklist files
blacklist_file = "../blacklists/blacklist_all.tsv"
mito_file = "../../../../data/AX4/genome_ver_2_7/gff3/chromosome_M.gff"
all_nc_file = "../../../../data/AX4/genome_ver_2_7/nc_genes.txt" # some overlap with IDs in 'blacklist_file', also includes 'splice transcripts'

# Paralog file
full_paralog_file = "../../../../data/AX4/genome_ver_2_7/biomart/genes_v52/d_disc_paralogs_biomart_ensembl_protist_ver_52.txt"

# Load paralog_data
full_paralog_data = CSV.read(full_paralog_file, DataFrame)
singleton_data = filter(row -> ismissing(row[15]), full_paralog_data)
select!(singleton_data, [1])
rename!(singleton_data, ["GeneID"])

# Load blacklist
blacklist = CSV.read(blacklist_file, DataFrame)

if filter_tes
    tes = filter(row -> row.type == "transposons" || row.type == "retro_transposons", blacklist).GeneID
    filter!(row -> row[1] ∉ tes, singleton_data)

end

if filter_ncRNAs
    ncRNAs = filter(row -> row.type == "transfer_rnas" || row.type == "ribosomal_rnas", blacklist).GeneID
    filter!(row -> row[1] ∉ ncRNAs, singleton_data)

end

if filter_all_nc
    all_nc_list = CSV.read(all_nc_file, DataFrame, header=false)[:,1]
    filter!(row -> row.GeneID ∉ all_nc_list, singleton_data)

end

if filter_plasmids
    plasmids = filter(row -> row.type == "plasmid", blacklist).GeneID
    filter!(row -> row[1] ∉ plasmids, singleton_data)

end

if filter_pseudogenes
    pseudogenes = filter(row -> row.type == "pseudogenes", blacklist).GeneID
    filter!(row -> row[1] ∉ pseudogenes, singleton_data)

end

if filter_ax4_dup_region
    ax4_dup_region = filter(row -> row.type == "ax4_dup", blacklist).GeneID
    filter!(row -> row[1] ∉ ax4_dup_region, singleton_data)

end

if filter_mito
    chrom_lengths_file = "../../../../data/AX4/genome_ver_2_7/chromosome_lengths.txt"
    include("./custom_lib/load_gff.jl")
    mito_ref = loadgff(mito_file, chrom_lengths_file; feature_type = "gene")
    mito_ids = [gene for gene in mito_ref.genes[1]]
    filter!(row -> row[1] ∉ mito_ids, singleton_data)

end

# Filter the two ensembl rna ids
filter!(row -> startswith(row.GeneID, "DDB_"), singleton_data)

# Save
CSV.write("./data/filtered/singleton_filt.tsv", singleton_data)