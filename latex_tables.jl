include("prelude.jl")

# LaTeX rows for the manuscript tables, built from the stats/ CSVs. Writes one file per table to
# ../../dicty_data/stats/latex/; paste the rows into dicty.tex / dicty_supplement.tex.

using Printf

const OUT = mkpath(joinpath(STATS_DIR, "latex"))
read_stats(name) = CSV.read(joinpath(STATS_DIR, "$name.csv"), DataFrame)

# p values: 3 decimals, scientific below 0.001, a bound where the value underflowed to 0
fmt_p(p) = ismissing(p) || isnan(p) ? "--" : p == 0 ? "\$<\$1e-300" : p < 0.001 ? @sprintf("%.2e", p) : @sprintf("%.3f", p)
fmt(x; d=3) = @sprintf("%.*f", d, round(x; digits=d) + 0.0)  # + 0.0 turns -0.0 into 0.0
bold(s, on) = on ? "\\textbf{$s}" : s
write_rows(name, rows) = open(io -> foreach(r -> println(io, r), rows), joinpath(OUT, "$name.tex"), "w")

const MARK_ROWS = [("H3K4me3", "H3K4me3"), ("H3K27ac", "H3K27ac"), ("H3K9me3", "H3K9me3"), ("ATAC", "Open Chromatin")]

# Table 2
function table2_rows(df)
    map(MARK_ROWS) do (mark, label)
        r = only(eachrow(filter(row -> row.Mark == mark, df)))
        "$label & $(fmt(r.r)) & $(bold(fmt_p(r.r_p_adj), r.r_p_adj < 0.05)) & $(fmt(r.eta2)) & $(bold(fmt_p(r.KW_p_adj), r.KW_p_adj < 0.05)) \\\\"
    end
end
write_rows("table2_expression", table2_rows(read_stats("enrichment_vs_expression_All")))
write_rows("table2_ds", table2_rows(read_stats("enrichment_vs_ds_All")))

# Table 3 and the per-stage regression table
const MLR_LABELS = [
    ("H3K4me3", "H3K4me3 (coverage)"), ("HasH3K4me3", "H3K4me3 (presence/absence)"),
    ("H3K27ac", "H3K27ac (coverage)"), ("HasH3K27ac", "H3K27ac (presence/absence)"),
    ("H3K9me3", "H3K9me3 (coverage)"), ("HasH3K9me3", "H3K9me3 (presence/absence)"),
    ("ATAC", "Open chromatin (coverage)"), ("HasATAC", "Open chromatin (presence/absence)"),
    ("ParalogDist", "Paralog distance"), ("SameChrom", "Same chromosome"),
    ("(Intercept)", "(Intercept)"), ("AvgExpr", "Expression (pair average)")]

function mlr_rows(name)
    df = read_stats(name)
    rows = String[]
    for (key, label) in MLR_LABELS
        i = findfirst(==(key), df.Predictor)
        isnothing(i) && continue  # no ATAC predictors for the vegetative stage
        sig = !isnan(df.pvalue_adj[i]) && df.pvalue_adj[i] < 0.05
        push!(rows, "$(bold(label, sig)) & $(bold(fmt(df.Coef[i]), sig)) & $(bold(fmt_p(df.pvalue_adj[i]), sig)) \\\\")
    end
    model = read_stats("$(name)_model")
    stat(s) = only(model.Value[model.Statistic .== s])
    push!(rows, "\\midrule", "\\textbf{R\\textsuperscript{2}} & \\textbf{$(fmt(stat("R2")))} & \\textbf{$(fmt_p(stat("F_p")))} \\\\")
    return rows
end
write_rows("table3", mlr_rows("mlr_All_ensembl"))
stage_rows = String[]
for (stage, title) in [("F", "Fruiting Body"), ("M", "Mound"), ("V", "Vegetative")]
    stage == "F" || push!(stage_rows, "\\midrule")
    push!(stage_rows, "\\multicolumn{3}{@{}l}{\\textbf{$title}} \\\\", "\\midrule")
    append!(stage_rows, mlr_rows("mlr_$(stage)_ensembl"))
end
write_rows("tableS_mlr_stage", stage_rows)

# VIFs
vif = read_stats("mlr_All_ensembl_vif")
write_rows("tableS_vif", ["$label & $(fmt(only(vif.VIF[vif.Predictor .== key]))) \\\\"
                          for (key, label) in vcat(MLR_LABELS[end:end], MLR_LABELS[1:10])])

# Table S1: Kruskal-Wallis and Dunn's tests by dS decile
kw = read_stats("kruskal_wallis_by_ds")
s1_rows = String[]
for (test, title) in [("h3k4me3_by_ds", "H3K4me3 vs.\\ dS (\\textit{D. discoideum})"),
                      ("atac_by_ds", "Open chromatin vs.\\ dS (\\textit{D. discoideum})"),
                      ("expression_by_ds", "Expression vs.\\ dS (\\textit{D. discoideum})"),
                      ("h3k9me3_by_ds_human", "H3K9me3 vs.\\ dS (human)")]
    k = only(eachrow(filter(row -> row.test == test, kw)))
    isempty(s1_rows) || push!(s1_rows, "\\midrule")
    push!(s1_rows, "\\multicolumn{4}{@{}p{0.95\\linewidth}}{\\textbf{Dunn's test, $title} (Kruskal--Wallis \$\\chi^2\$ = $(fmt(k.chi_squared; d=2)), \\textit{p} value $(fmt_p(k.p)))} \\\\", "\\midrule")
    for r in eachrow(read_stats("$(test)_dunn"))
        push!(s1_rows, "$(r.comparison) & $(fmt(r.Z)) & $(fmt_p(r.p)) & $(bold(fmt_p(r.p_adj), r.p_adj < 0.05)) \\\\")
    end
end
write_rows("tableS1", s1_rows)

# Principal component regression
pcr = read_stats("pca_regression_All")
loading_cols = ["AvgExpr", "H3K27ac", "H3K4me3", "H3K9me3", "ATAC", "ParalogDist",
                "HasH3K27ac", "HasH3K4me3", "HasH3K9me3", "HasATAC", "SameChrom"] .* "_loading"
write_rows("tableS_pcr", [join(vcat([r.Predictor, fmt(r.Coef), bold(fmt_p(r.pvalue_adj), r.pvalue_adj < 0.05),
                                     fmt(r.prop_var_explained)], [fmt(r[c]) for c in loading_cols]), " & ") * " \\\\"
                          for r in eachrow(pcr)])
