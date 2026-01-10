include("prelude.jl")

using Intervals


# Functions:
function contiguous_ones(signal::BitVector)
    in_stretch = false
    stretches = []
    start_ind = 0
    end_ind = 0
    for i in eachindex(signal)
        if signal[i]
            if !in_stretch
                in_stretch = true
                start_ind = i
            end
        else
            if in_stretch
                in_stretch = false
                end_ind = i - 1
                push!(stretches, (start_ind, end_ind))
            end
            in_stretch = false
        end
    end

    return stretches
end

function interval_df(signal::BitVector, name::String)
    stretches = contiguous_ones(signal)
    return DataFrame(Chromosome = [name for stretch in stretches], Start = [stretch[1] for stretch in stretches], End = [stretch[2] for stretch in stretches])
end

# TE files
te_id_file = "blacklists/ensembl_te_ids_with_predicted.tsv"
pred_te_coord_file = "../../dicty_data/rm_genomes/ensembl_52/full/results_v1/onecodetofindthemall/Dictyostelium_discoideum.dicty_2.7.dna.toplevel_aggregated_compat.csv"

# Genome data
gff_source = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"

# Peak files
chip_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_1_ensembl52/"

# Load genome data
ref_genome = loadgenome(gff_source, chrom_lengths_file)
te_ids = open(te_id_file) do file
    readlines(file)
end

# Load peak data
peak_files = readdir(chip_peak_file_dir, join=true) |> 
             filter(fn -> endswith(fn, ".narrowPeak")) |> 
             filter(fn -> contains(lowercase(fn), "k9me3")) |> 
             filter(fn -> !contains(fn, r"_S[AB]+_"))
peak_data = binpeaks(peak_files, chrom_lengths_file)
peak_dfs = [reduce(vcat, [interval_df(chrom.signal, chrom.name) for chrom in samp.chroms]) for samp in peak_data.samples]

# Get the known and predicted TE coordinates:
pred_te_df = CSV.read(pred_te_coord_file, DataFrame)
te_genes_list = [gene for gene in ref_genome.genes[2] if gene.id in te_ids]
te_gene_df = DataFrame([:Chromosome => [gene.scaffold.name for gene in te_genes_list], :Start => [gene.gene_start for gene in te_genes_list], :End => [gene.gene_end for gene in te_genes_list]])
combined_df = vcat(select(pred_te_df, [:Chromosome, :Start, :End]), te_gene_df)

# Get the union of the intervals for each chromosome for TEs:
te_interval_sets = Dict([chrom => IntervalSet([row.Start..row.End for row in eachrow(sub_df)]) for (sub_df, chrom) in zip(groupby(combined_df, :Chromosome), unique(combined_df.Chromosome))])
peak_interval_sets = [Dict([chrom => IntervalSet([row.Start..row.End for row in eachrow(sub_df)]) for (sub_df, chrom) in zip(groupby(peak_df, :Chromosome), unique(peak_df.Chromosome))]) for peak_df in peak_dfs]
all_keys = union(keys(te_interval_sets), union(keys(peak_interval_sets[1]), keys(peak_interval_sets[2]), keys(peak_interval_sets[3])))

# Get the union of the intervals for h3k9me3 peaks
peak_interval_sets_combined = Dict([k => union(union(haskey(peak_interval_sets[1], k) ? peak_interval_sets[1][k] : IntervalSet(), 
                                           haskey(peak_interval_sets[2], k) ? peak_interval_sets[2][k] : IntervalSet()), 
                                    haskey(peak_interval_sets[3], k) ? peak_interval_sets[3][k] : IntervalSet()) for k in all_keys])

# get the intersection of TE and peak regions
overlap_regions = Dict([k => haskey(te_interval_sets, k) && haskey(peak_interval_sets_combined, k) ? intersect(te_interval_sets[k], peak_interval_sets_combined[k]) : IntervalSet() for k in all_keys])
total_overlap = [Intervals.span.(overlap_regions[k].items) for k in all_keys]
total_overlap = sum(map(v -> isempty(v) ? 0 : sum(v), total_overlap))

total_bases = [(haskey(te_interval_sets, k) && !isempty(te_interval_sets[k].items) ? sum(Intervals.span.(te_interval_sets[k].items)) : 0, 
                haskey(peak_interval_sets_combined, k) && !isempty(peak_interval_sets_combined[k]) ? sum(Intervals.span.(peak_interval_sets_combined[k].items)) : 0) for k in all_keys]

te_bases = sum([pair[1] for pair in total_bases])
peak_bases = sum([pair[2] for pair in total_bases])

total_overlap / te_bases
total_overlap / peak_bases

