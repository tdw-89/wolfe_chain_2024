include("prelude.jl")

# Gene counts for Table 1. "Reference" counts are protein-coding genes in the Ensembl reference used in the
# analysis (D. discoideum: Ensembl Protists 52; human: Ensembl 99, excluding NMD candidates as in
# human_paralog_filt.jl), split by whether BioMart lists any paralog for the gene. "Analysed" counts are the
# singletons and paralog-pair members (dS <= 3) that enter the analyses.

data_dir = "../../dicty_data"
nonmissing_unique(v) = unique(collect(skipmissing(v)))

# D. discoideum
dicty_coding = Set(readlines(joinpath(data_dir, "AX4/genome_ver_2_7/ensembl_52/cds_ids.txt")))
dicty_biomart = CSV.read(joinpath(data_dir, "AX4/genome_ver_2_7/biomart/genes_v52/d_disc_paralogs_biomart_ensembl_protist_ver_52.txt"),
                         DataFrame, delim='\t')
dicty_with_paralog = Set(dicty_biomart[.!ismissing.(dicty_biomart[!, 15]), 1]) ∩ dicty_coding
dicty_pairs = filter(row -> row.dS <= 3, CSV.read(joinpath(data_dir, "filtered/paralog_filt.tsv"), DataFrame))
dicty_singletons = CSV.read(joinpath(data_dir, "filtered/singleton_filt.tsv"), DataFrame).GeneID
dicty_te_dist = CSV.read(joinpath(data_dir, "te_distance.csv"), DataFrame)

# Human
human_nmd = Set(readlines(joinpath(data_dir, "mammals/primates/h_sapiens/Ensembl_99/nmd_candidate_ids.txt")))
human_coding = setdiff(Set(readlines(joinpath(data_dir, "mammals/primates/h_sapiens/Ensembl_99/CDS_IDs.txt"))), human_nmd)
human_biomart = CSV.read(joinpath(data_dir, "mammals/primates/h_sapiens/Ensembl_99/paralog_info_GRCh38_ensembl_99.txt"), DataFrame)
human_with_paralog = Set(human_biomart[.!ismissing.(human_biomart[!, 11]), 1]) ∩ human_coding
human_pairs = filter(row -> row.dS <= 3, CSV.read(joinpath(data_dir, "filtered/human_paralog_info_filt.csv"), DataFrame))
human_singletons = CSV.read(joinpath(data_dir, "filtered/human_singletons_filt.csv"), DataFrame).GeneID
human_te_dist = CSV.read(joinpath(data_dir, "te_distance_human_TE.csv"), DataFrame)

counts = DataFrame(Species = String[], Category = String[], Reference = Int[], Analysed = Int[])
for (species, coding, with_paralog, pairs, singletons) in [
        ("D. discoideum", dicty_coding, dicty_with_paralog, dicty_pairs, dicty_singletons),
        ("H. sapiens", human_coding, human_with_paralog, human_pairs, human_singletons)]
    n_paralog_genes = length(nonmissing_unique(vcat(pairs.GeneID, pairs.ParalogID)))
    n_singletons = length(unique(singletons))
    push!(counts, (species, "Singletons", length(coding) - length(with_paralog), n_singletons))
    push!(counts, (species, "Paralogs", length(with_paralog), n_paralog_genes))
    push!(counts, (species, "Total", length(coding), n_singletons + n_paralog_genes))
end
write_stats("gene_counts", counts)

# Genes overlapping a TE, among the coding genes in the TE distance analysis (D. discoideum: non-TE,
# non-blacklisted coding genes; human: coding genes on the main chromosomes)
te_overlap = DataFrame(
    Species = ["D. discoideum", "H. sapiens"],
    GenesOverlappingTE = [count(isequal(0), skipmissing(df.Distance)) for df in (dicty_te_dist, human_te_dist)],
    GenesTested = [nrow(dicty_te_dist), nrow(human_te_dist)])
te_overlap.Percent = te_overlap.GenesOverlappingTE ./ te_overlap.GenesTested .* 100
write_stats("genes_overlapping_te", te_overlap)
