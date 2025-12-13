using Serialization

# Custom lib src:
include("./custom_lib/genomic_data.jl")

# Input bam directory
bam_dir = "D:/wang_et_al_new/ChIPseq_trimmed_bams_sorted_nodups/"

# Output directory
serial_dir = "../../dicty_data/julia_serialized/"

# Read in the bam file paths
input_bam_files = filter(fn -> contains(lowercase(fn), "input") && endswith(fn, ".bam"),readdir(bam_dir, join=true))

input_count_experiment = getallcountvectors(input_bam_files)
input_count_experiment = average_bam_replicate_groups(input_count_experiment)
GC.gc()

serialize(joinpath(serial_dir, "input_count_experiment.jls"), input_count_experiment)
