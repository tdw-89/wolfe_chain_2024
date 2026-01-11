include("prelude.jl")

using Serialization

using .EnrichmentUtils
using .RepeatUtils

reload_peak_data = true
te_type = "TE"

# Peak files
# chip_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_1_ensembl52/"
chip_peak_file_dir = "../../dicty_data/mammals/primates/h_sapiens/ENCODE_histone_mods/"

# Genome data
gff_source = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.99.gff3"
chrom_lengths_file = "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/chromosome_lengths.txt"
paralog_file = "../../dicty_data/filtered/human_paralog_info_filt.csv"
te_dist_file = "../../dicty_data/te_distance_human_$te_type.csv"

# Dosage sensitivity data:
significant_overlap_file = "../../dicty_data/lift_over_data_filtered.csv"
mismatch_file = "../../dicty_data/lift_over_data_mismatches.csv"

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
    serialize("../../dicty_data/julia_serialized/human_h3k9me3_exper.jls", peak_data)
else
    peak_data = deserialize("../../dicty_data/julia_serialized/human_h3k9me3_exper.jls")
end

addtogenes!(ref_genome, peak_data)

# Which genes have K9me3?
ids_with_k9me3 = []

for id in te_dist_df.GeneID

    gene = get(ref_genome, id)
    
    has_k9 = any([
        any(sum(getsiginrange(gene, GeneRange(REGION(), REGION(), 0, 0), ind)) > 0) for ind in eachindex(gene.samples)
    ])
    if has_k9
        push!(ids_with_k9me3, id)
    end
end

# TE distance distribution among filtered coding genes vs paralogs without H3K9me3:
paralog_ids = vcat(paralog_data.GeneID, paralog_data.ParalogID)
te_distances = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∉ paralog_ids && id ∉ ids_with_k9me3, te_dist_df.GeneID)]
te_distances_dups = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∈ paralog_ids && id ∉ ids_with_k9me3, te_dist_df.GeneID)]

# What proportion of non-paralog genes overlap a TE vs. paralog genes?
cont_table = [count(d -> d == 0, te_distances) count(d -> d == 0, te_distances_dups); 
              count(d -> d > 0, te_distances) count(d -> d > 0, te_distances_dups)]

# Are the two distributions significantly different? (A: No - LTR, <- 0.41264004228669854 -> adjust(Bonferroni) -> 1
#                                                     A: Yes - SINEs, <- 9.416765131506125e-5 -> adjust(Bonferroni) -> 0.000376670605260245
#                                                     A: Yes - LINEs, <- 0.003829015903780406 -> adjust(Bonferroni) -> 0.015316063615121623
#                                                     A: Yes - All TEs) <- 0.0017718660013885204 -> adjust(Bonferroni) -> 0.007087464005554081
pvalue(MannWhitneyUTest(te_distances, te_distances_dups))

# Are the proportions of genes overlapping a TE significantly different between paralogs and non-paralogs? (A: No - LTRs, -> 0.1764005443808332 -> adjust(Bonferroni) ->  0.7056021775233328
#                                                                                                           A: Yes - SINEs, -> 0.0003089060474037944 -> adjust(Bonferroni) -> 0.0012356241896151776
#                                                                                                           A: Yes - LINEs, -> 0.0018902814592321418 -> adjust(Bonferroni) -> 0.007561125836928567
#                                                                                                           A: Yes - All TEs) -> 0.0035883749357294805 -> adjust(Bonferroni) -> 0.014353499742917922
pvalue(FisherExactTest(cont_table[1, 1], cont_table[1, 2], cont_table[2, 1], cont_table[2, 2]))
adjust([0.1764005443808332, 0.0003089060474037944, 0.0018902814592321418, 0.0035883749357294805], BenjaminiHochberg())

# What is the distribution of TE distances among coding genes w/ vs. w/o H3K9me3?
# What about among just the filtered paralogs?
te_distances_k9me3 = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∈ ids_with_k9me3 && id ∉ paralog_ids, te_dist_df.GeneID)]
te_distances_no_k9me3 = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∉ ids_with_k9me3 && id ∉ paralog_ids, te_dist_df.GeneID)]

# Is the distribution of TE distances significantly different between genes with and without K9me3? (A: Yes - LTRs,
#                                                                                                    A: Yes - LINEs,
#                                                                                                    A: Yes - SINEs,
#                                                                                                    A: Yes - All TEs) -> 3.547289538993912e-60
plot([box(y=te_distances_k9me3, name="TE Distance w/ K9me3", marker_color="blue"), 
      box(y=te_distances, name="TE Distance w/o k9me3", marker_color="red")])
pvalue(MannWhitneyUTest(te_distances_k9me3, te_distances_no_k9me3))
cont_table = [count(d -> d == 0, te_distances_k9me3) count(d -> d == 0, te_distances_no_k9me3); 
              count(d -> d > 0, te_distances_k9me3) count(d -> d > 0, te_distances_no_k9me3)]

# Are the proportions of genes overlapping a TE significantly different between genes with and without K9me3? (A: Yes - LTRs,
#                                                                                                           A: Yes - SINEs,
#                                                                                                           A: Yes - LINEs,
#                                                                                                           A: Yes - All TEs) -> 3.9994659703972094e-66
pvalue(FisherExactTest(cont_table[1, 1], cont_table[1, 2], cont_table[2, 1], cont_table[2, 2]))

# Are the distributions of TE distances significantly different between paralogs with and without K9me3? (A: Yes - LTRs,
#                                                                                                         A: Yes - SINEs,
#                                                                                                         A: Yes - LINEs
#                                                                                                         A: Yes - All TEs) -> 1.1686502276252259e-45
te_distances_k9me3_dups = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∈ ids_with_k9me3 && id ∈ paralog_ids, te_dist_df.GeneID)]
te_distances_no_k9me3_dups = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> !(id ∈ ids_with_k9me3) && id ∈ paralog_ids, te_dist_df.GeneID)]
display(plot([box(y=te_distances_k9me3_dups, name="TE Distance w/ K9me3 (Paralogs)", marker_color="blue"), 
      box(y=te_distances_dups, name="TE Distance w/o k9me3 (Paralogs)", marker_color="red")],
      Layout(yaxis=attr(title="TE Distance"), plot_bgcolor="rgba(0,0,0,0)")))
pvalue(MannWhitneyUTest(te_distances_k9me3_dups, te_distances_no_k9me3_dups))

# Are the ratios of TE overlap significantly different between paralogs with and without K9me3? (A: Yes - LTRs,
#                                                                                                A: Yes - SINEs,
#                                                                                                A: Yes - LINEs,
#                                                                                                A: Yes - All TEs) -> 2.866036139619633e-48
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

# Are paralogs with H3K9me3 more likely to be triplosensitive? (A: No)
pvalue(MannWhitneyUTest(indiv_df.pTriplo[indiv_df.HasK9me3 .== true], 
                        indiv_df.pTriplo[indiv_df.HasK9me3 .== false]))

display(plot([
    box(y=indiv_df.pTriplo[indiv_df.TEDistance .<= 0], name="TE overlap", marker_color="blue"),
    box(y=indiv_df.pTriplo[indiv_df.TEDistance .> 0], name="No TE overlap", marker_color="red")
],
Layout(plot_bgcolor="rgba(0,0,0,0)",
yaxis=attr(gridcolor="lightgray",
           gridwidth=1.5))))

# Are genes with TE overlap more likely to be triplosensitive? (A: No - LTR -> 0.479435080559795,
#                                                               A: Yes - SINEs -> 1.2221709967990017e-8,
#                                                               A: Yes - LINEs -> 3.7744741554628688e-6,
#                                                               A: Yes - DNA TEs -> 1.0166622340827324e-9,
#                                                               A: Yes - All TEs -> 8.059948448185499e-6,
#                                                               )
pvalue(MannWhitneyUTest(indiv_df.pTriplo[indiv_df.TEDistance .<= 0], indiv_df.pTriplo[indiv_df.TEDistance .> 0]))

cont_table = [count(sens -> sens >= 0.5, indiv_df.pTriplo[indiv_df.TEDistance .<= 0]) count(sens -> sens < 0.5, indiv_df.pTriplo[indiv_df.TEDistance .<= 0]); 
              count(sens -> sens >= 0.5, indiv_df.pTriplo[indiv_df.TEDistance .> 0]) count(sens -> sens < 0.5, indiv_df.pTriplo[indiv_df.TEDistance .> 0])]

pvalue(FisherExactTest(cont_table[1, 1], cont_table[1, 2], cont_table[2, 1], cont_table[2, 2]))

# High-confidence sensitivity vs. TE distance vs. K9me3:
confidence_cutoff = 0.93

# Are high-confidence triplosensitive genes closer to TEs? (A: ? - LTR,
#                                                           A: Yes - SINEs,
#                                                           A: No - LINEs,
#                                                           A: Yes - All TEs) -> 0.03803689084188746
pvalue(MannWhitneyUTest(indiv_df.TEDistance[indiv_df.pTriplo .>= confidence_cutoff], indiv_df.TEDistance[indiv_df.pTriplo .< confidence_cutoff]))

# Are high-confidence triplosensitive genes more likely to have H3K9me3? (A: No)
cont_table1 = [count(sens -> sens >= confidence_cutoff, indiv_df.pTriplo[indiv_df.HasK9me3 .== true]) count(sens -> sens < confidence_cutoff, indiv_df.pTriplo[indiv_df.HasK9me3 .== true]); 
              count(sens -> sens >= confidence_cutoff, indiv_df.pTriplo[indiv_df.HasK9me3 .== false]) count(sens -> sens < confidence_cutoff, indiv_df.pTriplo[indiv_df.HasK9me3 .== false])]

sens_greater1 = (cont_table1[1, 1]/cont_table1[2,1]) > (cont_table1[1, 2]/cont_table1[2, 2])
p1 = pvalue(FisherExactTest(cont_table1[1, 1], cont_table1[1, 2], cont_table1[2, 1], cont_table1[2, 2]))

# Are high-confidence triplosensitive genes more likely to overlap a TE? (A: No - LTR,
#                                                                         A: Yes - SINEs,
#                                                                         A: Yes - LINEs,
#                                                                         A: Yes - All TEs) -> 0.02313507995779753
cont_table2 = [count(sens -> sens >= confidence_cutoff, indiv_df.pTriplo[indiv_df.TEDistance .<= 0]) count(sens -> sens < confidence_cutoff, indiv_df.pTriplo[indiv_df.TEDistance .<= 0]); 
              count(sens -> sens >= confidence_cutoff, indiv_df.pTriplo[indiv_df.TEDistance .> 0]) count(sens -> sens < confidence_cutoff, indiv_df.pTriplo[indiv_df.TEDistance .> 0])]

sens_greater2 = (cont_table2[1, 1]/cont_table2[2,1]) > (cont_table2[1, 2]/cont_table2[2, 2])
p2 = pvalue(FisherExactTest(cont_table2[1, 1], cont_table2[1, 2], cont_table2[2, 1], cont_table2[2, 2]))