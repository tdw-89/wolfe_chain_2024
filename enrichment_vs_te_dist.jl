include("prelude.jl")

using .EnrichmentUtils
using .RepeatUtils
using CategoricalArrays

# Mann–Whitney U tests of total TE density (upstream + gene body + downstream), one row each
test_results = DataFrame()

# Duplicate filtering parameter:
ds_threshold = 3

# Peak files
chip_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_1_ensembl52/"
atac_peak_file_dir = "../../dicty_data/wang_et_al/processed/run_2_ensembl52/"

# Genome data
gff_source = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/Dictyostelium_discoideum.dicty_2.7.52.gff3"
chrom_lengths_file = "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
paralog_file = "../../dicty_data/filtered/paralog_filt.tsv"
te_density_file = "../../dicty_data/te_density/te_density.csv"

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
                    yaxis=attr(title="TE Density",
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

# Load the TE density data:
te_density_df = CSV.read(te_density_file, DataFrame)
rename!(te_density_df, :Gene_ID => :GeneID)
te_density_df.TotalDensity = te_density_df.Upstream .+ te_density_df.Internal .+ te_density_df.Downstream

# Load the paralog data:
paralog_data = CSV.read(paralog_file, DataFrame)

# Filter the paralogs:
select!(paralog_data, ["GeneID", "ParalogID", "dS"])
filter!(row -> row["dS"] <= ds_threshold, paralog_data)

# Load peak data
peak_files = [readdir(chip_peak_file_dir, join=true); readdir(atac_peak_file_dir, join=true)]
peak_files = filter(fn -> endswith(fn, ".narrowPeak"), peak_files)
peak_files = filter(fn -> !contains(fn, r"_S[AB]+_"), peak_files)
peak_data = binpeaks(peak_files, chrom_lengths_file)
addtogenes!(ref_genome, peak_data)

# Which genes have H3K9me3?
ids_with_k9me3 = []
for id in te_density_df.GeneID
    gene = get(ref_genome, id)
    
    hasK9me3 = any([
        sum(getsiginrange(gene, GeneRange(TSS(), TES(), -500, 500), ind)) > 0 for ind in [3, 6, 9]
    ])
    if hasK9me3
        push!(ids_with_k9me3, id)
    end
end

# TE density distribution among filtered coding genes vs paralogs:
paralog_ids = vcat(paralog_data.GeneID, paralog_data.ParalogID)
low_ds_paralogs = paralog_data[levelcode.(cut(paralog_data.dS, 10)) .<= 2,:]
paralog_ids_low_ds = vcat(low_ds_paralogs.GeneID, low_ds_paralogs.ParalogID)
te_densities = te_density_df.TotalDensity[map(id -> id ∉ paralog_ids && id ∉ ids_with_k9me3, te_density_df.GeneID)]
te_densities_dups = te_density_df.TotalDensity[map(id -> id ∈ paralog_ids && id ∉ ids_with_k9me3, te_density_df.GeneID)]
fig = plot([box_plt_blue(te_densities, "TE Density"), box_plt_red(te_densities_dups, "TE Density (Paralogs)")],
      box_layout)
isinteractive() && display(fig)
savefig(fig, joinpath(@__DIR__, "../../dicty_data", "te_density_nonparalog_vs_paralog.html"))

push!(test_results, mwu_row("Non-paralogs vs paralogs (both without H3K9me3)", te_densities, te_densities_dups))

# What is the distribution of TE densities by H3K9me3 status?
te_densities_k9me3 = te_density_df.TotalDensity[map(id -> id ∈ ids_with_k9me3 && id ∉ paralog_ids, te_density_df.GeneID)]
te_densities_no_k9me3 = te_density_df.TotalDensity[map(id -> id ∉ ids_with_k9me3 && id ∉ paralog_ids, te_density_df.GeneID)]

fig = plot([box_plt_blue(te_densities_k9me3, "With H3K9me3"), box_plt_red(te_densities_no_k9me3, "Without H3K9me3")], box_layout)
isinteractive() && display(fig)
savefig(fig, joinpath(@__DIR__, "../../dicty_data", "te_density_by_k9me3_status.html"))
save_figure(plot([fig_box(te_densities_k9me3, "With H3K9me3", BLUE), fig_box(te_densities_no_k9me3, "Without H3K9me3", RED)],
                 fig_layout(yaxis=fig_axis(title="TE density"), xaxis=fig_axis(showgrid=false))), "figS06"; width=700, height=500)
push!(test_results, mwu_row("Non-paralogs with vs without H3K9me3", te_densities_k9me3, te_densities_no_k9me3))

te_densities_k9me3_dups = te_density_df.TotalDensity[map(id -> id ∈ ids_with_k9me3 && id ∈ paralog_ids, te_density_df.GeneID)]
te_densities_no_k9me3_dups = te_density_df.TotalDensity[map(id -> !(id ∈ ids_with_k9me3) && id ∈ paralog_ids, te_density_df.GeneID)]

fig = plot([box_plt_blue(te_densities_k9me3_dups, "With H3K9me3 (Paralogs)"), box_plt_red(te_densities_no_k9me3_dups, "Without H3K9me3 (Paralogs)")], box_layout)
isinteractive() && display(fig)
savefig(fig, joinpath(@__DIR__, "../../dicty_data", "te_density_paralogs_by_k9me3_status.html"))
save_figure(plot([fig_box(te_densities_k9me3_dups, "With H3K9me3", BLUE), fig_box(te_densities_no_k9me3_dups, "Without H3K9me3", RED)],
                 fig_layout(yaxis=fig_axis(title="TE density"), xaxis=fig_axis(showgrid=false))), "figS11"; width=700, height=500)
push!(test_results, mwu_row("Paralogs with vs without H3K9me3", te_densities_k9me3_dups, te_densities_no_k9me3_dups))

# Low dS paralogs
te_densities_low_ds = te_density_df.TotalDensity[map(id -> id ∈ paralog_ids_low_ds, te_density_df.GeneID)]
te_densities_low_ds_k9me3 = te_density_df.TotalDensity[map(id -> id ∈ paralog_ids_low_ds && id ∈ ids_with_k9me3, te_density_df.GeneID)]
te_densities_low_ds_no_k9me3 = te_density_df.TotalDensity[map(id -> id ∈ paralog_ids_low_ds && id ∉ ids_with_k9me3, te_density_df.GeneID)]

fig = plot([box_plt_blue(te_densities_low_ds_k9me3, "With H3K9me3"), box_plt_red(te_densities_low_ds_no_k9me3, "Without H3K9me3")], box_layout)
isinteractive() && display(fig)
savefig(fig, joinpath(@__DIR__, "../../dicty_data", "te_density_low_ds_paralogs_by_k9me3_status.html"))
push!(test_results, mwu_row("Low-dS paralogs with vs without H3K9me3", te_densities_low_ds_k9me3, te_densities_low_ds_no_k9me3))
test_results.p_adj = adjust(test_results.p, BenjaminiHochberg())
write_stats("te_density_tests", test_results)