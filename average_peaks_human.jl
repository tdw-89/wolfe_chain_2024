using CSV
using DataFrames
using Serialization

include("custom_lib/load_gff.jl")
include("custom_lib/genomic_data.jl")
include("custom_lib/misc_utils.jl")

reload_peak_data = false

chrom_lengths_file = "../../../../data/mammals/primates/h_sapiens/Ensembl_99/chromosome_lengths.txt"
peak_data_dir = "../../../../data/mammals/primates/h_sapiens/ENCODE_histone_mods/"

# Add the peak data
if reload_peak_data
    peak_files = reduce(vcat, [map(fn -> joinpath(root, fn), files) for (root, dir, files) in walkdir(chip_peak_file_dir)])
    peak_files = filter(fn -> endswith(fn, ".bed") || endswith(fn, ".bed.gz"), peak_files)
    peak_files = filter(fn -> contains(fn, "k9me3"), peak_files)
    peak_data = binpeaks(peak_files, chrom_lengths_file)
    serialize("./data/julia_serialized/human_h3k9me3_exper.jls", peak_data)
else
    peak_data = deserialize("./data/julia_serialized/human_h3k9me3_exper.jls")
end

peak_data_avg = average_peak_replicate_groups(peak_data)
bed_df = sample_to_bed_df(peak_data_avg.samples[1])
CSV.write("./data/human_h3k9me3_avg.bed", bed_df, delim='\t', header=false)