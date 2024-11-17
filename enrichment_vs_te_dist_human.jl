using CSV
using DataFrames
using PlotlyJS
using HypothesisTests
using MultipleTesting
using Serialization

# custom lib:
include("custom_lib/load_gff.jl")
include("custom_lib/enrichment_utils.jl")
include("custom_lib/te_utils.jl")

reload_peak_data = false
te_type = "TE"

# Peak files
# chip_peak_file_dir = "../../../../data/wang_et_al/processed/run_1_ensembl52/"
chip_peak_file_dir = "../../../../data/mammals/primates/h_sapiens/ENCODE_histone_mods/"

# Genome data
gff_source = "../../../../data/mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.99.gff3"
chrom_lengths_file = "../../../../data/mammals/primates/h_sapiens/Ensembl_99/chromosome_lengths.txt"
paralog_file = "./data/filtered/human_paralog_info_filt.csv"
te_dist_file = "./data/te_distance_human_$te_type.csv"

# Dosage sensitivity data: 
significant_overlap_file = "./data/lift_over_data_filtered.csv"
mismatch_file = "./data/lift_over_data_mismatches.csv"

# Load the genome data:
ref_genome = loadgenome(gff_source, chrom_lengths_file)

# Load the TE distance data:
te_dist_df = CSV.read(te_dist_file, DataFrame)

# Load the paralog data:
paralog_data = CSV.read(paralog_file, DataFrame)

# Filter the paralogs:
select!(paralog_data, ["GeneID", "ParalogID", "dS"])
filter!(row -> row.dS <= 3, paralog_data)

# Load peak data
if reload_peak_data
    peak_files = reduce(vcat, [map(fn -> joinpath(root, fn), files) for (root, dir, files) in walkdir(chip_peak_file_dir)])
    peak_files = filter(fn -> endswith(fn, ".bed") || endswith(fn, ".bed.gz"), peak_files)
    peak_files = filter(fn -> contains(fn, "k9me3"), peak_files)
    peak_data = binpeaks(peak_files, chrom_lengths_file)
    serialize("./data/julia_serialized/human_h3k9me3_exper.jls", peak_data)
else
    peak_data = deserialize("./data/julia_serialized/human_h3k9me3_exper.jls")
end

addtogenes!(ref_genome, peak_data)

# Which genes have K9me3?
ids_with_k9me3 = []

for id in te_dist_df.GeneID

    gene = get(ref_genome, id)
    
    has_k9 = any([
        any(getsiginrange(gene, GeneRange(REGION(), REGION(), 0, 0), ind)) for ind in eachindex(gene.samples)
    ])
    if has_k9
        push!(ids_with_k9me3, id)
    end
end

# TE distance distribution among filtered coding genes vs paralogs without H3K9me3:
paralog_ids = vcat(paralog_data.GeneID, paralog_data.ParalogID)
te_distances = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∉ paralog_ids && id ∉ ids_with_k9me3, te_dist_df.GeneID)]
te_distances_dups = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∈ paralog_ids && id ∉ ids_with_k9me3, te_dist_df.GeneID)]
plot([box(y=te_distances, name="TE Distance", marker_color="blue"), box(y=te_distances_dups, name="TE Distance (Paralogs)", marker_color="red")])

# What proportion of non-paralog genes overlap a TE vs. paralog genes?
cont_table = [count(d -> d == 0, te_distances) count(d -> d == 0, te_distances_dups); 
              count(d -> d > 0, te_distances) count(d -> d > 0, te_distances_dups)]

# Are the two distributions significantly different? (A: No - LTR,
#                                                     A: Yes - SINEs,
#                                                     A: Yes - LINEs,
#                                                     A: Yes - All TEs)
pvalue(MannWhitneyUTest(te_distances, te_distances_dups))

# Are the proportions of genes overlapping a TE significantly different between paralogs and non-paralogs? (A: No - LTRs
#                                                                                                           A: Yes - SINEs,
#                                                                                                           A: Yes - LINEs,
#                                                                                                           A: Yes - All TEs)
pvalue(FisherExactTest(cont_table[1, 1], cont_table[1, 2], cont_table[2, 1], cont_table[2, 2]))
pvalue(ChisqTest(cont_table))

# What is the distribution of TE distances among coding genes w/ vs. w/o H3K9me3?
# What about among just the filtered paralogs?
te_distances_k9me3 = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∈ ids_with_k9me3 && id ∉ paralog_ids, te_dist_df.GeneID)]
te_distances_no_k9me3 = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∉ ids_with_k9me3 && id ∉ paralog_ids, te_dist_df.GeneID)]

# Is the distribution of TE distances significantly different between genes with and without K9me3? (A: Yes - LTRs,
#                                                                                                    A: Yes - LINEs,
#                                                                                                    A: Yes - SINEs)
plot([box(y=te_distances_k9me3, name="TE Distance w/ K9me3", marker_color="blue"), box(y=te_distances, name="TE Distance w/o k9me3", marker_color="red")])
pvalue(MannWhitneyUTest(te_distances_k9me3, te_distances_no_k9me3))
cont_table = [count(d -> d == 0, te_distances_k9me3) count(d -> d == 0, te_distances_no_k9me3); 
              count(d -> d > 0, te_distances_k9me3) count(d -> d > 0, te_distances_no_k9me3)]

# Are the proportions of genes overlapping a TE significantly different between genes with and without K9me3? (A: Yes - LTRs,
#                                                                                                           A: Yes - SINEs,
#                                                                                                           A: Yes - LINEs,
#                                                                                                           A: Yes - All TEs)
pvalue(FisherExactTest(cont_table[1, 1], cont_table[1, 2], cont_table[2, 1], cont_table[2, 2]))

# Are the distributions of TE distances significantly different between paralogs with and without K9me3? (A: Yes - LTRs,
#                                                                                                         A: Yes - SINEs,
#                                                                                                         A: Yes - LINEs)
te_distances_k9me3_dups = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∈ ids_with_k9me3 && id ∈ paralog_ids, te_dist_df.GeneID)]
te_distances_no_k9me3_dups = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> !(id ∈ ids_with_k9me3) && id ∈ paralog_ids, te_dist_df.GeneID)]
display(plot([box(y=te_distances_k9me3_dups, name="TE Distance w/ K9me3 (Paralogs)", marker_color="blue"), 
      box(y=te_distances_dups, name="TE Distance w/o k9me3 (Paralogs)", marker_color="red")],
      Layout(yaxis=attr(title="TE Distance"), plot_bgcolor="rgba(0,0,0,0)")))
pvalue(MannWhitneyUTest(te_distances_k9me3_dups, te_distances_no_k9me3_dups))

# Are the ratios of TE overlap significantly different between paralogs with and without K9me3? (A: Yes - LTRs,
#                                                                                                A: Yes - SINEs,
#                                                                                                A: Yes - LINEs,
#                                                                                                A: Yes - All TEs)
plot(bar(x=["Paralogs With H3K9me3", "Paralogs Without H3K9me3"], 
         y=[count(d -> d == 0, te_distances_k9me3_dups)/length(te_distances_k9me3_dups), 
            count(d -> d == 0, te_distances_no_k9me3_dups)/length(te_distances_no_k9me3_dups)]),
            Layout(yaxis=attr(title="% Overlapping a TE", range=[0, 1])))
cont_table = [count(d -> d == 0, te_distances_k9me3_dups) count(d -> d == 0, te_distances_no_k9me3_dups); 
              count(d -> d > 0, te_distances_k9me3_dups) count(d -> d > 0, te_distances_no_k9me3_dups)]

pvalue(FisherExactTest(cont_table[1, 1], 
                       cont_table[1, 2], 
                       cont_table[2, 1], 
                       cont_table[2, 2]))

# Triplosensitivity vs. TE distance vs. K9me3:
significant_overlap_df = CSV.read(significant_overlap_file, DataFrame)
select!(significant_overlap_df, ["gene_id_GRCh38", "gene_name_liftover_GRCh37", "pTriplo_liftover_GRCh37"])
rename!(significant_overlap_df, Dict(:gene_id_GRCh38 => :GeneID, :gene_name_liftover_GRCh37 => :Gene, :pTriplo_liftover_GRCh37 => :pTriplo))

# Individual gene sensitivity vs. TE distance vs. K9me3:
indiv_df = DataFrame(
    :GeneID => vcat(paralog_data.GeneID, paralog_data.ParalogID)
)
indiv_df = innerjoin(indiv_df, significant_overlap_df, on=:GeneID)

insertcols!(indiv_df, :TEDistance => 0.0, :HasK9me3 => false)
indiv_df.TEDistance = [only(te_dist_df.Distance[te_dist_df.GeneID .== gene_id]) for gene_id in indiv_df.GeneID]
indiv_df.HasK9me3 = [gene_id ∈ ids_with_k9me3 for gene_id in indiv_df.GeneID]

display(plot([
    box(y=indiv_df.pTriplo[indiv_df.HasK9me3 .== true], name="With K9me3", marker_color="blue"),
    box(y=indiv_df.pTriplo[indiv_df.HasK9me3 .== false], name="Without K9me3", marker_color="red")
],
Layout(plot_bgcolor="rgba(0,0,0,0)",
yaxis=attr(gridcolor="lightgray",
           gridwidth=1.5))
))

# Are genes with H3K9me3 more likely to be triplosensitive? (A: No)
pvalue(MannWhitneyUTest(indiv_df.pTriplo[indiv_df.HasK9me3 .== true], indiv_df.pTriplo[indiv_df.HasK9me3 .== false]))

display(plot([
    box(y=indiv_df.pTriplo[indiv_df.TEDistance .<= 0], name="TE overlap", marker_color="blue"),
    box(y=indiv_df.pTriplo[indiv_df.TEDistance .> 0], name="No TE overlap", marker_color="red")
],
Layout(plot_bgcolor="rgba(0,0,0,0)",
yaxis=attr(gridcolor="lightgray",
           gridwidth=1.5))))

# Are genes with TE overlap more likely to be triplosensitive? (A: No - LTR,
#                                                               A: Yes - SINEs,
#                                                               A: Yes - LINEs,
#                                                               A: Yes - All TEs)
pvalue(MannWhitneyUTest(indiv_df.pTriplo[indiv_df.TEDistance .<= 0], indiv_df.pTriplo[indiv_df.TEDistance .> 0]))

cont_table = [count(sens -> sens >= 0.5, indiv_df.pTriplo[indiv_df.TEDistance .<= 0]) count(sens -> sens < 0.5, indiv_df.pTriplo[indiv_df.TEDistance .<= 0]); 
              count(sens -> sens >= 0.5, indiv_df.pTriplo[indiv_df.TEDistance .> 0]) count(sens -> sens < 0.5, indiv_df.pTriplo[indiv_df.TEDistance .> 0])]

pvalue(FisherExactTest(cont_table[1, 1], cont_table[1, 2], cont_table[2, 1], cont_table[2, 2]))

# High-confidence sensitivity vs. TE distance vs. K9me3:
confidence_cutoff = 0.93

# Are high-confidence triplosensitive genes closer to TEs? (A: ? - LTR,
#                                                           A: Yes - SINEs,
#                                                           A: No - LINEs,
#                                                           A: ? - All TEs)
pvalue(MannWhitneyUTest(indiv_df.TEDistance[indiv_df.pTriplo .>= confidence_cutoff], indiv_df.TEDistance[indiv_df.pTriplo .< confidence_cutoff]))

# Are high-confidence triplosensitive genes more likely to have H3K9me3? (A: No)
cont_table1 = [count(sens -> sens >= confidence_cutoff, indiv_df.pTriplo[indiv_df.HasK9me3 .== true]) count(sens -> sens < confidence_cutoff, indiv_df.pTriplo[indiv_df.HasK9me3 .== true]); 
              count(sens -> sens >= confidence_cutoff, indiv_df.pTriplo[indiv_df.HasK9me3 .== false]) count(sens -> sens < confidence_cutoff, indiv_df.pTriplo[indiv_df.HasK9me3 .== false])]

sens_greater1 = (cont_table1[1, 1]/cont_table1[2,1]) > (cont_table1[1, 2]/cont_table1[2, 2])
p1 = pvalue(FisherExactTest(cont_table1[1, 1], cont_table1[1, 2], cont_table1[2, 1], cont_table1[2, 2]))

# Are high-confidence triplosensitive genes more likely to overlap a TE? (A: No - LTR,
#                                                                         A: Yes - SINEs,
#                                                                         A: Yes - LINEs,
#                                                                         A: Yes - All TEs)
cont_table2 = [count(sens -> sens >= confidence_cutoff, indiv_df.pTriplo[indiv_df.TEDistance .<= 0]) count(sens -> sens < confidence_cutoff, indiv_df.pTriplo[indiv_df.TEDistance .<= 0]); 
              count(sens -> sens >= confidence_cutoff, indiv_df.pTriplo[indiv_df.TEDistance .> 0]) count(sens -> sens < confidence_cutoff, indiv_df.pTriplo[indiv_df.TEDistance .> 0])]

sens_greater2 = (cont_table2[1, 1]/cont_table2[2,1]) > (cont_table2[1, 2]/cont_table2[2, 2])
p2 = pvalue(FisherExactTest(cont_table2[1, 1], cont_table2[1, 2], cont_table2[2, 1], cont_table2[2, 2]))

counts = [count(sens -> sens >= confidence_cutoff, indiv_df.pTriplo[indiv_df.HasK9me3 .== true .&& indiv_df.TEDistance .<= 0]),
          count(sens -> sens < confidence_cutoff, indiv_df.pTriplo[indiv_df.HasK9me3 .== true .&& indiv_df.TEDistance .<= 0]),
          count(sens -> sens >= confidence_cutoff, indiv_df.pTriplo[indiv_df.HasK9me3 .== false .&& indiv_df.TEDistance .<= 0]), 
          count(sens -> sens < confidence_cutoff, indiv_df.pTriplo[indiv_df.HasK9me3 .== false .&& indiv_df.TEDistance .<= 0]),
          count(sens -> sens >= confidence_cutoff, indiv_df.pTriplo[indiv_df.HasK9me3 .== true .&& indiv_df.TEDistance .> 0]),
          count(sens -> sens < confidence_cutoff, indiv_df.pTriplo[indiv_df.HasK9me3 .== true .&& indiv_df.TEDistance .> 0]),
          count(sens -> sens >= confidence_cutoff, indiv_df.pTriplo[indiv_df.HasK9me3 .== false .&& indiv_df.TEDistance .> 0]),
          count(sens -> sens < confidence_cutoff, indiv_df.pTriplo[indiv_df.HasK9me3 .== false .&& indiv_df.TEDistance .> 0])]

display(plot([
    bar(x=["Hi-Confidence With K9me3 & TE overlap", 
            "Low-Confidence With K9me3 & TE overlap",
            "Hi-Confidence Without K9me3 & TE overlap",
            "Low-Confidence Without K9me3 & TE overlap",
            "Hi-Confidence With K9me3 & No TE overlap",
            "Low-Confidence With K9me3 & No TE overlap",
            "Hi-Confidence Without K9me3 & No TE overlap",
            "Low-Confidence Without K9me3 & No TE overlap"],
        y=counts)
]))

sens_greater3 = (cont_table3[1, 1]/cont_table3[2,1]) > (cont_table3[1, 2]/cont_table3[2, 2])
p3 = pvalue(FisherExactTest(cont_table3[1, 1], cont_table3[1, 2], cont_table3[2, 1], cont_table3[2, 2]))

df = DataFrame(
    :Var => ["K9me3", "TE Distance", "K9me3 & TE Distance"],
    :SensGreater => [sens_greater1, sens_greater2, sens_greater3],
    :PValue => [p1, p2, p3]
)


# OLD:
#=
# Interaction count vs. TE distance vs. K9me3:
interact_count_df = CSV.read(prot_interact_count_file, DataFrame)
insertcols!(paralog_data, :InteractCountGene => 0, :InteractCountParalog => 0)

for i in 1:nrow(paralog_data)
    gene_id = paralog_data.GeneID[i]
    paralog_id = paralog_data.ParalogID[i]
    
    gene_interact_count = interact_count_df.count[interact_count_df.gene .== gene_id]
    paralog_interact_count = interact_count_df.count[interact_count_df.gene .== paralog_id]
    
    if !isempty(gene_interact_count)
        paralog_data.InteractCountGene[i] = only(gene_interact_count)
    else
        paralog_data.InteractCountGene[i] = -1
    end
    
    if !isempty(paralog_interact_count)
        paralog_data.InteractCountParalog[i] = only(paralog_interact_count)
    else
        paralog_data.InteractCountParalog[i] = -1
    end
end

both_counts = filter(row -> row.InteractCountGene != -1 && row.InteractCountParalog != -1, paralog_data)
insertcols!(both_counts, :AvgCount => 0.0, :TEDistance => 0.0, :HasK9me3 => false)

# Caclulate the average interaction count:
for i in 1:nrow(both_counts)
    both_counts.AvgCount[i] = mean([both_counts.InteractCountGene[i], both_counts.InteractCountParalog[i]])
    both_counts.TEDistance[i] = only(te_dist_df.Distance[te_dist_df.GeneID .== both_counts.GeneID[i]])
    both_counts.HasK9me3[i] = both_counts.GeneID[i] ∈ ids_with_k9me3
end

display(plot([
    box(y=both_counts.AvgCount[both_counts.HasK9me3 .== true], name="With K9me3", marker_color="blue"),
    box(y=both_counts.AvgCount[both_counts.HasK9me3 .== false], name="Without K9me3", marker_color="red")
]))

pvalue(MannWhitneyUTest(both_counts.AvgCount[both_counts.HasK9me3 .== true], both_counts.AvgCount[both_counts.HasK9me3 .== false]))

display(plot([
    box(y=both_counts.AvgCount[both_counts.TEDistance .<= 0], name="TE overlap", marker_color="blue"),
    box(y=both_counts.AvgCount[both_counts.TEDistance .> 0], name="No TE overlap", marker_color="red")
]))

pvalue(MannWhitneyUTest(both_counts.AvgCount[both_counts.TEDistance .<= 0], both_counts.AvgCount[both_counts.TEDistance .> 0]))
=#
