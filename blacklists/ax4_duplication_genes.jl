using CSV
using DataFrames
using FastaIO
using Intervals

cd("blacklists")

ensembl_cds = "../../../../../data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.cds.all.fa"

const ax4_dup_interval_1 = Interval(2263132, 3015703)
const ax4_dup_interval_2 = Interval(3016083, 3768654)

# Get all Ensembl CDS sequences:
ensembl_reader = FastaReader(ensembl_cds)
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

insertcols!(ensembl_cds_df, :interval => [Interval(ensembl_cds_df.Start[i], ensembl_cds_df.End[i]) for i in 1:nrow(ensembl_cds_df)])

ensembl_cds_df_chrom2 = filter(row -> row.Chrom == "2", ensembl_cds_df)

# Get the genes that overlap with either of the two intervals:
overlapping_genes = filter(row -> Intervals.overlaps(row.interval, ax4_dup_interval_1) || Intervals.overlaps(row.interval, ax4_dup_interval_2), ensembl_cds_df_chrom2)

# Write the genes to a file:
CSV.write("ax4_duplication_genes.tsv", DataFrame(overlapping_genes))