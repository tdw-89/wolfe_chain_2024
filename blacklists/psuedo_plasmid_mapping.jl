using CSV
using DataFrames
using FastaIO
using Intervals
using BioSequences
using BioAlignments

cd("./blacklists")

##############
# FUNCTIONS: #
##############
affinegap = AffineGapScoreModel(
    match=5,
    mismatch=-4,
    gap_open=-5,
    gap_extend=-3
)
const comp_dict = Dict(['A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', 'N' => 'N'])
function reverse_complement(seq::String)
    return join([comp_dict[base] for base in reverse(seq)])
end
function deg_overlap(interval1, interval2)

    if Intervals.overlaps(interval1, interval2)
        # What % of interval1 is contained within interval2?
        int1_in_int2 = (span(Intervals.intersect(interval1, interval2)) / span(interval1)) * 100
        # What % of interval2 is contained within interval1?
        int2_in_int1 = (span(Intervals.intersect(interval1, interval2)) / span(interval2)) * 100

        return int1_in_int2, int2_in_int1
    else
        return 0, 0
    end
end
function align_seqs(seq1, seq2, score_model)
    seq1 = LongDNA{4}(seq1)
    seq2 = LongDNA{4}(seq2)
    return pairalign(GlobalAlignment(), seq1, seq2, score_model)
end

# Genome data:
gene_info_file = "../../dicty_data/AX4/gene_information.txt"
ensembl_cds = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.cds.all.fa"
ensembl_ncrna = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.ncrna.fa"
dictybase_cds = "../../dicty_data/AX4/genome_ver_2_7/fastas/dicty_primary_cds.fa"
dictybase_nc_file = "../../dicty_data/AX4/genome_ver_2_7/fastas/dicty_noncoding.fa"
chrom_names_file = "../../dicty_data/AX4/genome_ver_2_7/chromosome_names.csv"

# Load the gene info data
gene_info = CSV.read(gene_info_file, DataFrame)
rename!(gene_info, ["GeneID", "GeneName", "Synonyms", "GeneProducts"])

# Load the chromosome synonyms:
chrom_names = CSV.read(chrom_names_file, DataFrame)
chrom_name_dict = Dict(chrom_names[!,1] .=> chrom_names[!,2])

# Get the plasmids:
plasmids = filter(row -> contains(row.GeneName, "Ddp") ||
                         (!ismissing(row.Synonyms) && contains(lowercase(row.Synonyms), "plasmid")) ||
                         (!ismissing(row.GeneProducts) && contains(lowercase(row.GeneProducts), "plasmid")), gene_info)

# Get the pseudogenes:
pseudogenes = filter(row -> contains(row.GeneName, "_ps") || (!ismissing(row.GeneProducts) && contains(lowercase(row.GeneProducts), "pseudogene")), gene_info)

# Load the CDS data:
ensembl_reader = FastaReader(ensembl_cds)

# Get all Ensembl CDS sequences:
descs = []
seqs = []
for (desc, seq) in ensembl_reader
    # if split(strip(match(r"gene:[A-Z0-9a-z_]+ ", desc).match), ":")[2] in missing_te_ids
    push!(descs, desc)
    push!(seqs, seq)
end

ensembl_cds_df = DataFrame(GeneID = String[], Chrom = String[], Start = Int[], End = Int[], Length = Int[], Strand = Char[], Sequence = String[])
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
    push!(ensembl_cds_df, [id, chrom, start_pos, end_pos, seq_length, strand, seq])
end

# Get all Ensembl non-coding RNA sequences:
ensembl_nc_reader = FastaReader(ensembl_ncrna)
descs = []
seqs = []
for (desc, seq) in ensembl_nc_reader
    push!(descs, desc)
    push!(seqs, seq)
end

ensembl_nc_df = DataFrame(GeneID = String[], Chrom = String[], Start = Int[], End = Int[], Length = Int[], Strand = Char[], Sequence = String[])
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
    push!(ensembl_nc_df, [id, chrom, start_pos, end_pos, seq_length, strand, seq])
end

# Get all dictybase CDS sequences:
dictybase_cds_reader = FastaReader(dictybase_cds)
descs = []
seqs = []
for (desc, seq) in dictybase_cds_reader
    push!(descs, desc)
    push!(seqs, seq)
end

dictybase_cds_df = DataFrame(GeneID = String[], Chrom = String[], Start = Int[], End = Int[], Length = Int[], Sequence = String[])
for (i, desc) in enumerate(descs)
    desc_split = split(desc, "|")
    id = strip(desc_split[2])
    if !isnothing(match(r"chromosome:.*", desc))
        pos_str = strip(split(match(r"chromosome:.*", desc).match,":")[2])
        pos_split = split(pos_str)
        chrom = strip(pos_split[1])
        start_pos = parse(Int, pos_split[3])
        end_pos = parse(Int, pos_split[5])
    else
        chrom = "NA"
        start_pos = -1
        end_pos = -1
    end
    seq_length = end_pos - start_pos + 1
    seq = seqs[i]
    push!(dictybase_cds_df, [id, chrom, start_pos, end_pos, seq_length > 1 ? seq_length : -1, seq])
end

# Get all dictybase non-coding sequences:
dictybase_nc_reader = FastaReader(dictybase_nc_file)
descs = []
seqs = []
for (desc, seq) in dictybase_nc_reader
    push!(descs, desc)
    push!(seqs, seq)
end

dictybase_nc_df = DataFrame(GeneID = String[], Chrom = String[], Start = Int[], End = Int[], Length = Int[], Sequence = String[])
for (i, desc) in enumerate(descs)
    desc_split = split(desc, "|")
    id = strip(desc_split[2])
    pos_str = strip(split(match(r"chromosome:.*", desc).match,":")[2])
    pos_split = split(pos_str)
    chrom = strip(pos_split[1])
    start_pos = parse(Int, pos_split[3])
    end_pos = parse(Int, pos_split[5])
    seq_length = end_pos - start_pos + 1
    seq = seqs[i]
    push!(dictybase_nc_df, [id, chrom, start_pos, end_pos, seq_length, seq])
end

# Get the positions of the pseudogenes and plasmids:
pseudogenes_df = filter(row -> row.GeneID in pseudogenes.GeneID, dictybase_nc_df)
@assert nrow(pseudogenes_df) == nrow(pseudogenes)
plasmids_df = filter(row -> row.GeneID in plasmids.GeneID, dictybase_cds_df)
@assert nrow(plasmids_df) == nrow(plasmids)

# Which plasmids are missing from the Ensembl CDS data?
missing_plasmids = filter(row -> row.GeneID ∉ ensembl_cds_df.GeneID, plasmids_df)

# Since all are missing, do any Ensembl CDS sequences match plasmid sequences?
findall(map(seq -> seq in plasmids_df.Sequence, ensembl_cds_df.Sequence))

# Which pseudogenes are missing from the Ensembl CDS data?
matching_pseudogenes = innerjoin(pseudogenes_df, ensembl_cds_df, on = :GeneID, makeunique=true)
missing_pseudogenes = filter(row -> row.GeneID ∉ ensembl_cds_df.GeneID, pseudogenes_df)

#####################
# THOSE THAT MATCH: #
#####################

    # Which are on the same chromosome with the same coordinates?

same_chrom_same_coords = filter(row -> row.Chrom == row.Chrom_1 && row.Start == row.Start_1 && row.End == row.End_1, matching_pseudogenes)

    # Which are on the same chromosome with different coordinates?

same_chrom_diff_coords = filter(row -> row.Chrom == row.Chrom_1 && (row.Start != row.Start_1 || row.End != row.End_1), matching_pseudogenes)

    # Which are on different chromosomes with the same sequence?

diff_chrom_same_seq = filter(row -> row.Chrom != row.Chrom_1 && row.Sequence == row.Sequence_1, matching_pseudogenes)

    # Which are on different chromosomes with different sequences?

diff_chrom_diff_seq = filter(row -> row.Chrom != row.Chrom_1 && row.Sequence != row.Sequence_1, matching_pseudogenes)

# For those with the same chromosome but different coordinates, how much do they overlap?
insertcols!(same_chrom_diff_coords, 7, :DbInterval => [Interval(row.Start, row.End) for row in eachrow(same_chrom_diff_coords)])
insertcols!(same_chrom_diff_coords, :EnsemblInterval => [Interval(row.Start_1, row.End_1) for row in eachrow(same_chrom_diff_coords)])
insertcols!(same_chrom_diff_coords, 8 , :OverlapEnsPercent => zeros(nrow(same_chrom_diff_coords)))
insertcols!(same_chrom_diff_coords, :OverlapDbPercent => zeros(nrow(same_chrom_diff_coords)))
for i in 1:nrow(same_chrom_diff_coords)
    db_interval = same_chrom_diff_coords.DbInterval[i]
    ens_interval = same_chrom_diff_coords.EnsemblInterval[i]
    db_overlap_ens, ens_overlap_db = deg_overlap(db_interval, ens_interval)
    same_chrom_diff_coords.OverlapEnsPercent[i] = db_overlap_ens
    same_chrom_diff_coords.OverlapDbPercent[i] = ens_overlap_db
end

insertcols!(same_chrom_diff_coords, :AvgOverlap => (same_chrom_diff_coords.OverlapEnsPercent .+ same_chrom_diff_coords.OverlapDbPercent) ./ 2)

# Which Ensembl genes are totally overlapped by a Db pseudogene?
same_chrom_diff_coords_100 = filter(row -> row.OverlapDbPercent == 100, same_chrom_diff_coords)

# Which Ensembl genes overlap a Db pseudogene by at least 10%?
same_chrom_diff_coords_10 = filter(row -> row.OverlapDbPercent >= 10, same_chrom_diff_coords)

# Are there any pseudogenes that don't overlap with Ensembl CDS?
zero_overlap = filter(row -> row.OverlapEnsPercent == 0 || row.OverlapDbPercent == 0, same_chrom_diff_coords)

# Are there any Ensembl CDS genes that overlap with Db pseudogenes?
insertcols!(ensembl_cds_df, :interval => [Interval(row.Start, row.End) for row in eachrow(ensembl_cds_df)])
insertcols!(pseudogenes_df, :interval => [Interval(row.Start, row.End) for row in eachrow(pseudogenes_df)])
non_pseudo_cds_ens = filter(row -> row.GeneID ∉ pseudogenes_df.GeneID, ensembl_cds_df)
overlapping_ids_ens = []
overlapped_ids_pseudo = []
for i in 1:nrow(non_pseudo_cds_ens)
    cds_ens_interval = non_pseudo_cds_ens.interval[i]
    cds_ens_chrom = non_pseudo_cds_ens.Chrom[i]
    temp_db_pseudo_df = filter(row -> row.Chrom == cds_ens_chrom, pseudogenes_df)
    for k in 1:nrow(temp_db_pseudo_df)
        if Intervals.overlaps(cds_ens_interval, temp_db_pseudo_df.interval[k])
            push!(overlapping_ids_ens, non_pseudo_cds_ens.GeneID[i])
            push!(overlapped_ids_pseudo, temp_db_pseudo_df.GeneID[k])
        end
    end
end
overlapping_cds_ens = map(id -> findfirst(id .== non_pseudo_cds_ens.GeneID), overlapping_ids_ens)
overlapping_cds_ens = non_pseudo_cds_ens[overlapping_cds_ens, :]
overlapped_pseudo = map(id -> findfirst(id .== pseudogenes_df.GeneID), overlapped_ids_pseudo)
overlapped_pseudo = pseudogenes_df[overlapped_pseudo, :]
overlapping_df = hcat(overlapping_cds_ens, overlapped_pseudo, makeunique=true)

# What % and how many nucleotides of each overlapping Ensembl CDS sequence is contained within the overlapping pseudogene?
insertcols!(overlapping_df, 9, :OverlapPercentEns => zeros(nrow(overlapping_df)))
insertcols!(overlapping_df, 10, :OverlappingNucEns => zeros(Int, nrow(overlapping_df)))

for i in 1:nrow(overlapping_df)
    ens_overlap = intersect(overlapping_df.interval[i], overlapping_df.interval_1[i])
    overlapping_nuc = span(ens_overlap)
    overlap_percent = (overlapping_nuc / span(overlapping_df.interval[i])) * 100
    overlapping_df.OverlapPercentEns[i] = overlap_percent
    overlapping_df.OverlappingNucEns[i] = overlapping_nuc
end

# Which overlapping Ensembl CDS sequences have an overlap % == 100?
overlapping_df_100 = filter(row -> row.OverlapPercentEns == 100, overlapping_df)

# Which overlapping Ensembl CDS sequences overlap with a Db pseudogene by at least 10%?
overlapping_df_10 = filter(row -> row.OverlapPercentEns >= 10, overlapping_df)

# Which overlapping Ensembl CDS sequences have an overlap % less than 100%, start with a start codon, end with a stop codon, and have a sequence length divisible by 3?
overlapping_df_lessthan100 = filter(row -> row.OverlapPercentEns < 100, overlapping_df)
overlapping_df_lessthan100 = filter(row -> startswith(row.Sequence, "ATG") && (endswith(row.Sequence, "TAA") || endswith(row.Sequence, "TAG") || endswith(row.Sequence, "TGA")) && length(row.Sequence) % 3 == 0, overlapping_df_lessthan100)

###########################
# THOSE THAT DON'T MATCH: #
###########################

# Are any of the missing IDs in the Ensembl ncRNA fasta? A: NO
missing_ncrna = filter(row -> row.GeneID in missing_pseudogenes.GeneID, ensembl_nc_df)

# Write the pseudogenes to a file:

                        # These match perfectly, so keep
pseudogene_list = vcat(same_chrom_same_coords.GeneID,
                        # These Ens genes are partially overlapped by the same Db pseudogene
                        same_chrom_diff_coords_10.GeneID,
                        # This gene has the same sequence but is on a non-numbered chrom in both:
                        diff_chrom_same_seq.GeneID, 
                        # The Ens and Db sequences differ only by small overhangs
                        diff_chrom_diff_seq.GeneID,
                        overlapping_df.GeneID)


CSV.write("ensembl_pseudogenes.tsv", DataFrame(GeneID=pseudogene_list))


# OLD ↓↓↓↓

#=

# Do any Ensembl coding genes overlap any missing dictybase pseudogenes that are on numbered chromosomes?
insertcols!(missing_pseudogenes, :interval => [Interval(row.Start, row.End) for row in eachrow(missing_pseudogenes)])
overlapping_ids_ens = []
for i in 1:nrow(ensembl_cds_df)

    ensembl_id = ensembl_cds_df.GeneID[i]
    ensembl_chrom = ensembl_cds_df.Chrom[i]

    if ensembl_chrom in ["$i" for i in 1:6]
        ensembl_interval = Interval(ensembl_cds_df.Start[i], ensembl_cds_df.End[i])
        for j in 1:nrow(missing_pseudogenes)
            if missing_pseudogenes.Chrom[j] == ensembl_chrom
                if Intervals.overlaps(ensembl_interval, missing_pseudogenes.interval[j])
                    println("Ensembl gene $ensembl_id on chromosome $ensembl_chrom overlaps missing pseudogene $(missing_pseudogenes.GeneID[j])")
                    push!(overlapping_ids_ens, ensembl_id)
                end
            end
        end
    end
end

overlapping_missing_pseudogenes = filter(row -> row.GeneID in overlapping_ids_ens, ensembl_cds_df)

# Among the pseudogene IDs that are present, do they have the same sequence?
matching_join_df = innerjoin(pseudogenes_df, ensembl_cds_df, on = :GeneID, makeunique=true)
no_pseudo_match_ens_cds = filter(row -> row.GeneID ∉ pseudogenes_df.GeneID, ensembl_cds_df)

# Do any coding genes that don't match dictybase pseudogene IDs overlap any dictybase pseudogenes that are on numbered chromosomes?
insertcols!(no_pseudo_match_ens_cds, :interval => [Interval(row.Start, row.End) for row in eachrow(no_pseudo_match_ens_cds)])
insertcols!(pseudogenes_df, :interval => [Interval(row.Start, row.End) for row in eachrow(pseudogenes_df)])
overlaps_pseudo = []
is_overlapped = []
for i in 1:nrow(no_pseudo_match_ens_cds)

    ensembl_id = no_pseudo_match_ens_cds.GeneID[i]
    ensembl_chrom = no_pseudo_match_ens_cds.Chrom[i]
    ensembl_interval = no_pseudo_match_ens_cds.interval[i]

    if ensembl_chrom in ["$i" for i in 1:6]
        
        temp_pseudo_df = filter(row -> row.Chrom == ensembl_chrom, pseudogenes_df)
        for j in 1:nrow(temp_pseudo_df)
            if Intervals.overlaps(ensembl_interval, temp_pseudo_df.interval[j])
                # println("Ensembl gene $ensembl_id on chromosome $ensembl_chrom overlaps pseudogene $(pseudogenes_df.GeneID[j])")
                push!(overlaps_pseudo, ensembl_id)
                push!(is_overlapped, temp_pseudo_df.GeneID[j])
            end
        end
    end
end

overlap_df_ens = [findfirst(id .== no_pseudo_match_ens_cds.GeneID) for id in overlaps_pseudo]
overlap_df_ens = no_pseudo_match_ens_cds[overlap_df_ens, :]
overlap_df_pseudo = [findfirst(id .== pseudogenes_df.GeneID) for id in is_overlapped]
overlap_df_pseudo = pseudogenes_df[overlap_df_pseudo, :]
overlap_df_combined = hcat(overlap_df_ens, overlap_df_pseudo, makeunique=true)



same_coords = filter(row -> row.Start == row.Start_1 && row.End == row.End_1 && row.Chrom == row.Chrom_1, matching_join_df)
same_coords_same_seq = filter(row -> row.Sequence == row.Sequence_1, same_coords)
same_coords_diff_seq = filter(row -> row.Sequence != row.Sequence_1, same_coords)


# TESTING #
ens_genomic = raw"C:\Users\eccwo\Dropbox\UML\Grad_school\research\data\AX4\genome_ver_2_7\ensembl_52\Dictyostelium_discoideum.dicty_2.7.dna.toplevel.fa.gz"
ens_genomic_reader = FastaReader(ens_genomic)
db_genomic = raw"C:\Users\eccwo\Dropbox\UML\Grad_school\research\data\AX4\genome_ver_2_7\fastas\dicty_chromosomal.fa"
db_genomic_reader = FastaReader(db_genomic)

# Ensembl genomic data:
ens_descs = []
ens_seqs = []
for (desc, seq) in ens_genomic_reader
    push!(ens_descs, desc)
    push!(ens_seqs, seq)
end
ens_genomic_df = DataFrame(Desc = ens_descs, Sequence = ens_seqs)

# Dictybase genomic data:
db_descs = []
db_seqs = []
for (desc, seq) in db_genomic_reader
    push!(db_descs, desc)
    push!(db_seqs, seq)
end
db_genomic_df = DataFrame(Desc = db_descs, Sequence = db_seqs)

db_genomic_df[5,2] == ens_genomic_df[1,2] # Chromosome 1
db_genomic_df[6,2] == ens_genomic_df[2,2] # Chromosome 2
db_genomic_df[7,2] == ens_genomic_df[3,2] # Chromosome 3
db_genomic_df[8,2] == ens_genomic_df[4,2] # Chromosome 4
db_genomic_df[9,2] == ens_genomic_df[5,2] # Chromosome 5
db_genomic_df[10,2] == ens_genomic_df[6,2] # Chromosome 6

reverse_complement(ens_genomic_df[5,2][35600:40649])

=#