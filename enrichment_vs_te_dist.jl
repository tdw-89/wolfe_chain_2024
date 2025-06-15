using CSV
using DataFrames
using PlotlyJS
using HypothesisTests
using MultipleTesting

# custom lib:
include("custom_lib/load_gff.jl")
include("custom_lib/enrichment_utils.jl")
include("custom_lib/te_utils.jl")

# Duplicate filtering parameter:
id_method = "dS" #
id_threshold = 3

# Peak files
# chip_peak_file_dir = "../../../../data/wang_et_al/processed/run_1_ensembl52/"
chip_peak_file_dir = "../../../../data/wang_et_al/processed/run_1_ensembl52/"
atac_peak_file_dir = "../../../../data/wang_et_al/processed/run_2_ensembl52/"

# Genome data
gff_source = "../../../../data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../../../data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
paralog_file = "./data/filtered/paralog_filt.tsv"
final_gene_list = "./data/filtered/final_gene_list.txt"
te_dist_file = "./data/te_distance.csv"
# te_dist_file = "./data/repeat_distance.csv"

font_family = "Times New Roman"
box_plt_blue(data, name) = box(y=data, 
                        name=name, 
                        marker_color="navy", 
                        line=attr(color="navy"), 
                        fillcolor="#596fff")
box_plt_red(data, name) = box(y=data, 
                        name=name, 
                        marker_color="red", 
                        line=attr(color="red"), 
                        fillcolor="#fc8b8b")
box_layout = Layout(titlefont=attr(family=font_family, size=40),
                    plot_bgcolor="white",
                    yaxis=attr(title="log(TE Distance)",
                                titlefont=attr(family=font_family, size=30),
                                tickfont=attr(family=font_family, size=30),
                                gridcolor="lightgray",
                                gridwidth=1.5), 
                    xaxis=attr(titlefont=attr(family=font_family, size=30),
                                tickfont=attr(family=font_family, size=30),
                                gridcolor="lightgray",
                                gridwidth=1.5), 
                    showlegend=false)

# Load the genome data:
ref_genome = loadgenome(gff_source, chrom_lengths_file)

# Load the TE distance data:
te_dist_df = CSV.read(te_dist_file, DataFrame)
# te_dist_df.Distance = log10.(te_dist_df.Distance .+ 1)

# Load the paralog data:
paralog_data = CSV.read(paralog_file, DataFrame)

# Load final filtered gene list:
final_gene_list = open(final_gene_list) do file
    readlines(file)
end

# Filter the paralogs:
select!(paralog_data, ["GeneID", "ParalogID", id_method])
filter!(row -> row[id_method] <= id_threshold, paralog_data)

# Load peak data
peak_files = [readdir(chip_peak_file_dir, join=true); readdir(atac_peak_file_dir, join=true)]
peak_files = filter(fn -> endswith(fn, ".narrowPeak"), peak_files)
peak_files = filter(fn -> !contains(fn, r"_S[AB]+_"), peak_files)
peak_data = binpeaks(peak_files, chrom_lengths_file)
significant_regions = CSV.read("./data/sig_regions.csv", DataFrame)
addtogenes!(ref_genome, peak_data)

# Which genes have H3K9me3?
ids_with_k9me3 = []
k9me3_start_offset = only(significant_regions[significant_regions.Mark .== "K9me3",:Start])
k9me3_end_offset = only(significant_regions[significant_regions.Mark .== "K9me3",:End])

for id in te_dist_df.GeneID
    gene = get(ref_genome, id)
    
    hasK9me3 = any([
        sum(getsiginrange(gene, GeneRange(TSS(), TES(), k9me3_start_offset, k9me3_end_offset), ind)) > 0 for ind in [3, 6, 9]
    ])
    if hasK9me3
        push!(ids_with_k9me3, id)
    end
end

# TE distance distribution among filtered coding genes vs paralogs:
paralog_ids = vcat(paralog_data.GeneID, paralog_data.ParalogID)
low_ds_paralogs = paralog_data[levelcode.(cut(paralog_data.dS, 10)) .<= 2,:]
paralog_ids_low_ds = vcat(low_ds_paralogs.GeneID, low_ds_paralogs.ParalogID)
te_distances = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∉ paralog_ids && id ∉ ids_with_k9me3, te_dist_df.GeneID)]
te_distances_dups = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&& map(id -> id ∈ paralog_ids && id ∉ ids_with_k9me3, te_dist_df.GeneID)]
plot([box_plt_blue(log10.(te_distances), "TE Distance"), box_plt_red(log10.(te_distances_dups), "TE Distance (Paralogs)")],
      box_layout)

# What proportion of non-paralog genes are within 10kb of a TE vs. paralog genes?
cont_table = [count(d -> d <= 10_000, te_distances) count(d -> d <= 10_000, te_distances_dups); 
              count(d -> d > 10_000, te_distances) count(d -> d > 10_000, te_distances_dups)]

# Are the means of the two distributions significantly different? (A: No)
pvalue(FisherExactTest(cont_table[1, 1], cont_table[1, 2], cont_table[2, 1], cont_table[2, 2]))
pvalue(MannWhitneyUTest(te_distances, te_distances_dups))

# What is the distribution of TE distances among paralogs?
te_distances_k9me3 = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&&  map(id -> id ∈ ids_with_k9me3 && id ∉ paralog_ids, te_dist_df.GeneID)]
te_distances_no_k9me3 = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&&  map(id -> id ∉ ids_with_k9me3 && id ∉ paralog_ids, te_dist_df.GeneID)]

plot([box_plt_blue(log.(te_distances_k9me3), "With H3K9me3"), box_plt_red(log.(te_distances), "Without H3K9me3")], box_layout)
pvalue(MannWhitneyUTest(log.(te_distances_k9me3), log.(te_distances_no_k9me3)))

te_distances_k9me3_dups = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&&  map(id -> id ∈ ids_with_k9me3 && id ∈ paralog_ids, te_dist_df.GeneID)]
te_distances_no_k9me3_dups = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&&  map(id -> !(id ∈ ids_with_k9me3) && id ∈ paralog_ids, te_dist_df.GeneID)]

plot([box_plt_blue(log.(te_distances_k9me3_dups), "With H3K9me3 (Paralogs)"), box_plt_red(log.(te_distances_dups), "Without H3K9me3 (Paralogs)")], box_layout)
pvalue(MannWhitneyUTest(log.(te_distances_k9me3_dups), log.(te_distances_no_k9me3_dups)))

# Are the distributions of TE distances significantly different between paralogs with and without H3K9me3?
plot(bar(x=["Paralogs With H3K9me3", "Paralogs Without H3K9me3"], 
         y=[count(d -> d <= 10_000, te_distances_k9me3_dups)/length(te_distances_k9me3_dups), 
            count(d -> d <= 10_000, te_distances_no_k9me3_dups)/length(te_distances_no_k9me3_dups)]),
            Layout(yaxis=attr(title="% 10kb or less from TE", range=[0, 1])))
cont_table = [count(d -> d <= 10_000, te_distances_k9me3_dups) count(d -> d <= 10_000, te_distances_no_k9me3_dups); 
              count(d -> d > 10_000, te_distances_k9me3_dups) count(d -> d > 10_000, te_distances_no_k9me3_dups)]

pvalue(FisherExactTest(cont_table[1, 1], cont_table[1, 2], cont_table[2, 1], cont_table[2, 2]))

# Low dS paralogs
te_distances_low_ds = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&&  map(id -> id ∈ paralog_ids_low_ds, te_dist_df.GeneID)]
te_distances_low_ds_k9me3 = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&&  map(id -> id ∈ paralog_ids_low_ds && id ∈ ids_with_k9me3, te_dist_df.GeneID)]
te_distances_low_ds_no_k9me3 = te_dist_df.Distance[te_dist_df.Distance .!= Inf .&&  map(id -> id ∈ paralog_ids_low_ds && id ∉ ids_with_k9me3, te_dist_df.GeneID)]

plot([box_plt_blue(log.(te_distances_low_ds_k9me3), "With H3K9me3"), box_plt_red(log.(te_distances_low_ds_no_k9me3), "Without H3K9me3")], box_layout)