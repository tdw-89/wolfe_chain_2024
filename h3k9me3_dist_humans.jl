include("prelude.jl")

using Serialization

include("custom_lib/genomic_data.jl")

exon_coverage_file = "../../dicty_data/human_h3k9me3_exon_enrichment.csv"
chrom_lengths_file = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/chromosome_lengths.txt"
te_coverage_file = "../../dicty_data/h3k9me3_coverage_all.csv"
peak_data_file = "../../dicty_data/julia_serialized/human_h3k9me3_exper.jls"

# Load the h3k9me3 peak data
peak_data = deserialize(peak_data_file)

# Load the exon coverage data
exon_coverage = CSV.read(exon_coverage_file, DataFrame)

# Load the TE coverage data
repeat_coverage = CSV.read(te_coverage_file, DataFrame)

# Average the peak data:
peak_data = average_peak_replicate_groups(peak_data)

total_h3k9me3_bases = 
mapreduce(
    chr -> count(i -> i > 0, chr.signal),
    +,
    peak_data.samples[1].chroms
)

# Calculate the total number of bases in numbered chromosomes:
chrom_len_df = CSV.read(chrom_lengths_file, DataFrame, header=false)
filter!(row -> row[1] in [["$i" for i in 1:22]..., "X", "Y"], chrom_len_df)
total_chrom_bases = sum(chrom_len_df[!,2])

# Fraction of bases covered by h3k9me3 in at least one sample:
fraction_h3k9me3 = total_h3k9me3_bases / total_chrom_bases

# Fraction of H3K9me3 peaks present in at least one sample that cover an exon:
exon_coverage_n = mapreduce(
    pair -> round(pair[1] * pair[2]),
    +,
    zip(exon_coverage.Length, exon_coverage.OverlapPerc)
)

fraction_exon_coverage = exon_coverage_n / total_h3k9me3_bases

# Fraction of H3K9me3 peaks present in at least one sample that cover a TE:
TE_types = ["LINE", "SINE", "DNA", "Retroposon", "LTR", "RC"]
te_coverage = filter(row -> row.Type in TE_types, repeat_coverage)
te_coverage_n = mapreduce(
    pair -> round(pair[1] * pair[2]),
    +,
    zip(te_coverage.TotalBases, te_coverage.Coverage)
)

fraction_te_coverage = te_coverage_n / total_h3k9me3_bases

# Fraction of H3K9me3 peaks present in at least one sample that cover a non-TE repeat:
non_te_coverage = filter(row -> row.Type ∉ [TE_types..., "CDS", "Unknown"], repeat_coverage)

non_te_coverage_n = mapreduce(
    pair -> round(pair[1] * pair[2]),
    +,
    zip(non_te_coverage.TotalBases, non_te_coverage.Coverage)
)

fraction_non_te_coverage = non_te_coverage_n / total_h3k9me3_bases

# Fraction of H3K9me3 peaks present in at least one sample that cover a CDS:
cds_coverage_n = only(repeat_coverage.Coverage[repeat_coverage.Type .== "CDS"]) * only(repeat_coverage.TotalBases[repeat_coverage.Type .== "CDS"])
fraction_cds_coverage = cds_coverage_n / total_h3k9me3_bases

# Fraction of genic peaks covering exons:
genic_exon_coverage = exon_coverage_n / only(repeat_coverage.TotalBases[repeat_coverage.Type .== "CDS"])
