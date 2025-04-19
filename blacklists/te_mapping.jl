using CSV
using DataFrames
using FastaIO
using Intervals

cd("./blacklists")

##############
# FUNCTIONS: #
##############

const comp_dict = Dict(['A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', 'N' => 'N'])
function reverse_complement(seq::String)
    return join([comp_dict[base] for base in reverse(seq)])
end

# Custom lib src:
# include("./custom_lib/load_gff.jl")
# include("./custom_lib/genomic_data.jl")
# include("./custom_lib/misc_utils.jl")

# Genomic data
dictybase_cds = "../../../../../data/AX4/genome_ver_2_7/fastas/dicty_primary_cds.fa"
ensembl_cds = "../../../../../data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.cds.all.fa"
ensembl_paralogs_filt = "../data/filtered/paralog_filt.tsv"

# CDS Readers:
db_reader = FastaReader(dictybase_cds)
ensembl_reader = FastaReader(ensembl_cds)

##############################
# GET THE RELEVANT CDS INFO: #
##############################

# Get the dictybase annotated TE info:
descs = []
seqs = []
for (desc, seq) in db_reader
    if contains(desc, r"_[R]?TE")
        push!(descs, desc)
        push!(seqs, seq)
    end
end

db_te_df = DataFrame(GeneID = String[], Chrom = String[], Start = Int[], End = Int[], Length = Int[], Sequence = String[])
for (i, desc) in enumerate(descs)
    desc_split = split(desc, "|")
    id = strip(desc_split[2])
    chrom = match(r"chromosome: [0-9A-Z]+", desc).match
    chrom = strip(split(chrom)[2])
    pos = match(r"position .*", desc).match
    start_pos = parse(Int, split(pos)[2])
    end_pos = parse(Int, split(pos)[4])
    seq_length = end_pos - start_pos + 1
    seq = seqs[i]
    push!(db_te_df, [id, chrom, start_pos, end_pos, seq_length, seq])
end

# Get all Ensembl CDS sequences:
descs = []
seqs = []
for (desc, seq) in ensembl_reader
    # if split(strip(match(r"gene:[A-Z0-9a-z_]+ ", desc).match), ":")[2] in missing_te_ids
    push!(descs, desc)
    push!(seqs, seq)
end

ensembl_cds_df = DataFrame(GeneID = String[], Chrom = String[], Start = Int[], End = Int[], Length = Int[], Sequence = String[], Strand = Char[])
for (i, desc) in enumerate(descs)
    id = split(strip(match(r"gene:[A-Z0-9a-z_]+ ", desc).match), ":")[2]
    pos = strip(match(r"dicty_2.7:[A-Z0-9a-z:-]+ ", desc).match)
    chrom = split(pos, ":")[2]
    start_pos = parse(Int, split(pos, ":")[3])
    end_pos = parse(Int, split(pos, ":")[4])
    seq_length = end_pos - start_pos + 1
    strand = split(pos, ":")[5]
    strand = strand == "-1" ? '-' : strand == "1" ? '+' : '?'
    seq = seqs[i]
    push!(ensembl_cds_df, [id, chrom, start_pos, end_pos, seq_length, seq, strand])
end

###########################
# MAP THE TE ANNOTATIONS: #
###########################

# Get all Ensembl genes annotated as TEs by dictybase:
ensembl_te_df = filter(row -> row.GeneID in db_te_df.GeneID, ensembl_cds_df)

# Which TE IDs are missing from the ensembl data?
missing_tes = filter(row -> row.GeneID in setdiff(db_te_df.GeneID, ensembl_te_df.GeneID), db_te_df)

# Combine the Ensembl CDS and dictybase TE dataframes:
combined_te_df = innerjoin(db_te_df, ensembl_te_df, on=:GeneID, makeunique=true)

    # Which TEs are on the same chromosome with the same coordinates?

same_chrom_same_coord = filter(row -> row.Chrom == row.Chrom_1 && row.Start == row.Start_1 && row.End == row.End_1, combined_te_df)

    # Which TEs are on the same chromosome with different coordinates?

same_chrom_diff_coord = filter(row -> row.Chrom == row.Chrom_1 && (row.Start != row.Start_1 || row.End != row.End_1), combined_te_df)

    # Which TEs are on different chromosomes with the same sequence?

same_seq_diff_chrom = filter(row -> row.Chrom != row.Chrom_1 && row.Sequence == row.Sequence_1, combined_te_df)

    # Which TEs are on different chromosomes with different sequences?

diff_seq_diff_chrom = filter(row -> row.Chrom != row.Chrom_1 && row.Sequence != row.Sequence_1, combined_te_df)

# How much do TEs that are on the same chrom but have different coordinates overlap?
insertcols!(same_chrom_diff_coord, "Overlap %" => zeros(nrow(same_chrom_diff_coord)))
for row in eachrow(same_chrom_diff_coord)
    min_start = min(row.Start, row.Start_1)
    max_end = max(row.End, row.End_1)
    interval_len = max_end - min_start + 1
    len_1 = row.End - row.Start + 1
    len_2 = row.End_1 - row.Start_1 + 1
    overlap = (len_1 + len_2 - interval_len) / interval_len
    row."Overlap %" = overlap
end

# Are there any other coding genes that overlap with these dictybase (↑) TEs in the Ensembl reference (checked manually via Jbrowse)?
# Yes:
#   - DDB_G0271004 overlaps with the region spanned by the union of DDB_G0271006's coordinates in the both references
#   - DDB_G0277517 overlaps with the region spanned by the union of DDB_G0277519's coordinates in the both references
#   - DDB_G0267368 overlaps with the region spanned by the union of DDB_G0267366's coordinates in the both references
#   NOTE:
#   - DDB_G0271440 is just much shorter in the Ensembl reference, but the region spanned by this gene in dictybase does not contain any other coding genes in the Ensembl reference

# Which TEs are on different chromosomes and are on chromosomes 1-6 in one or both genomes?
diff_num_chroms = filter(row -> row.Chrom in ["$i" for i in 1:6] || row.Chrom_1 in ["$i" for i in 1:6], diff_seq_diff_chrom)

# Is there a coding gene that overlaps with this (↑) TE's dictybase coordinates in the Ensembl reference (checked manually via Jbrowse)? 
# No, there are no coding genes that overlap with this TE's dictybase coordinates in the Ensembl reference

# Do any Ensembl coding genes overlap any missing dictybase TEs? # 
insertcols!(missing_tes, :interval => [Interval(row.Start, row.End) for row in eachrow(missing_tes)])
overlapping_ids_ens = []
overlapping_ids_te = []
for i in 1:nrow(ensembl_cds_df)

    ensembl_id = ensembl_cds_df.GeneID[i]
    ensembl_chrom = ensembl_cds_df.Chrom[i]

    if ensembl_chrom in ["$i" for i in 1:6]
        ensembl_interval = Interval(ensembl_cds_df.Start[i], ensembl_cds_df.End[i])
        for j in 1:nrow(missing_tes)
            if missing_tes.Chrom[j] == ensembl_chrom
                if Intervals.overlaps(ensembl_interval, missing_tes.interval[j])
                    println("Ensembl gene $ensembl_id on chromosome $ensembl_chrom overlaps missing TE $(missing_tes.GeneID[j])")
                    push!(overlapping_ids_ens, ensembl_id)
                    push!(overlapping_ids_te, missing_tes.GeneID[j])
                end
            end
        end
    end
end

# Create a list of genes that overlap with missing or mis-mapped TEs from the dictybase reference:
overlapping_genes = unique(vcat(overlapping_ids_ens, 
                                [
                                    "DDB_G0271004",
                                    "DDB_G0277517",
                                    "DDB_G0267368",
                                ]))

# Create the full list of Ensembl TEs:
                           # Keep those that have the same positions on the same chromosome
te_ids_ensembl_full = vcat(same_chrom_same_coord.GeneID, 
                           # Keep those on non-numbered chromosomes with the same sequence (the one on a numbered chromosome in db didn't overlap any Ensembl gene or have the same sequence as the Ensembl gene that matched that ID) 
                           filter(row -> row.Chrom ∉ ["$i" for i in 1:6] && row.Chrom_1 ∉ ["$i" for i in 1:6], same_seq_diff_chrom).GeneID,
                           # Keep those on numbered chromosomes with different positions (they all overlap)
                           same_chrom_diff_coord.GeneID,
                           # Keep those on numbered chromosomes that overlap the 'same_chrom_diff_coord' TEs
                           overlapping_genes)

CSV.write("ensembl_te_ids_full.tsv", DataFrame(GeneID=te_ids_ensembl_full))

####################
# PARALOG OVERLAP: #
####################

# Load the paralog data:
paralog_data = CSV.read(ensembl_paralogs_filt, DataFrame)

# Do any of the paralog IDs match any of the Ensembl TE IDs?
paralog_ids = vcat(paralog_data.GeneID, paralog_data.ParalogID)
intersect(te_ids_ensembl_full, paralog_ids)

#= 
NOTES:
    - There were 10 TE IDs that were not present in the Ensembl reference
        - While three of these TEs overlapped with 2 Ensembl coding genes, neither of the Ensembl genes were in the semi-filtered (%ID > 30 && rbh) paralog data
    - Among the 510 IDs that were present, there were 452 that were on the same chromosome in both references
        - Same chrom: 
            - Among these, 4 had different start and end positions: DDB_G0271006, DDB_G0277519, DDB_G0271440, DDB_G0267366
        - Different chrom:
            - Among the 58 that were on different chromosomes, only one ID mapped to a numbered chromosome in one reference
=#