include("prelude.jl")

# Composite manuscript figures, built from the plot traces serialized by enrichment_vs_expression.jl,
# enrichment_vs_ds.jl and enrichment_vs_ds_human.jl and from the stats/ tables. Writes PDFs to FIG_DIR.
# Single-panel figures are exported by the scripts that compute them.

using CategoricalArrays
using Printf
using Serialization

ser_dir = "../../dicty_data/julia_serialized"
data_dir = "../../dicty_data"

# Serialized mark order is H3K27ac, H3K4me3, H3K9me3, ATAC; panels follow the captions
const MARK_NAMES = ["H3K27ac", "H3K4me3", "H3K9me3", "Open chromatin (ATAC)"]
const PANEL_ORDER = [2, 1, 3, 4]

# Region tick positions within the TSS (-500..0), gene-body (percentile) and TES (0..500) heatmaps
const TSS_TICKS = -400:100:-100
const BODY_TICKS = 20:20:80
const TES_TICKS = 100:100:400

copy_trace(t) = GenericTrace(deepcopy(t.fields))

# Heatmaps are embedded in the PDF as images, which viewers smooth when enlarged; repeating each decile row
# keeps the row boundaries sharp. Returns the trace on a numeric y axis and the tick positions and labels.
const ROW_UPSAMPLE = 20
function upsample_rows!(t)
    labels = collect(t[:y])
    m = Matrix{Float64}(t[:z])
    # repeat along whichever dimension holds the deciles; the trace's own orientation is kept
    t[:z] = size(m, 1) == length(labels) ? repeat(m, inner=(ROW_UPSAMPLE, 1)) : repeat(m, inner=(1, ROW_UPSAMPLE))
    t[:y] = collect(1:length(labels) * ROW_UPSAMPLE) .- 0.5
    return (collect(((1:length(labels)) .- 0.5) .* ROW_UPSAMPLE), labels)
end

load_traces(prefix, kind, stage) =
    [deserialize(joinpath(ser_dir, "$(prefix)$(region)_$(kind)$stage.jls"))[1]
     for region in ("tss_enrich_plots", "body_enrich_plots", "tes_enrich_plots", "bar_plots")]

# Dashed vertical guide lines in a heatmap
guide_lines(xref, yref, vals) = [attr(type="line", x0=v, x1=v, y0=0, y1=1, xref=xref, yref="$yref domain",
                                      line=attr(color="white", width=1, dash="dot")) for v in vals]

# Heatmaps of z-scored peak frequency (TSS | gene body | TES) with the box plot of the per-gene average
# beside them, one row per mark
function enrichment_figure(traces, marks, var_label, box_color)
    tss, body, tes, boxes = traces
    n = length(marks)
    gap = 0.05
    h = (1 - (n - 1) * gap) / n
    cols = [[0.0, 0.13], [0.135, 0.395], [0.40, 0.53], [0.67, 1.0]]
    data = GenericTrace[]
    layout = fig_layout(margin=attr(l=130, r=20, t=40, b=110))
    shapes, annotations = [], []
    tick_pos = Dict{Int, Any}()
    for (r, m) in enumerate(marks)
        k = 4 * (r - 1)
        y_dom = [1 - (r - 1) * (h + gap) - h, 1 - (r - 1) * (h + gap)]
        bottom = r == n
        for (j, t) in enumerate((tss[m], body[m], tes[m], boxes[m]))
            t2 = copy_trace(t)
            t2[:xaxis], t2[:yaxis] = "x$(k + j)", "y$(k + j)"
            if j <= 3
                row_ticks = upsample_rows!(t2)
                j == 1 && (tick_pos[r] = row_ticks)
                t2[:colorscale] = HEATMAP_COLORS
                t2[:showscale] = r == 1 && j == 3
                t2[:colorbar] = attr(x=0.54, xanchor="left", y=0.5, yanchor="middle", len=0.6, thickness=12,
                                     title=attr(text="z", side="right"), outlinewidth=0)
            else
                t2[:marker_color] = box_color
                t2[:line] = attr(color=box_color, width=1.5)
                t2[:marker] = attr(size=3, color=box_color)
            end
            push!(data, t2)
        end
        layout["xaxis$(k + 1)"] = attr(domain=cols[1], anchor="y$(k + 1)", range=[-500, 0], tickvals=collect(TSS_TICKS),
                                       ticktext=string.(TSS_TICKS), showticklabels=bottom, ticks="")
        layout["xaxis$(k + 2)"] = attr(domain=cols[2], anchor="y$(k + 2)", tickvals=collect(BODY_TICKS),
                                       ticktext=string.(BODY_TICKS) .* "%", showticklabels=bottom, ticks="")
        layout["xaxis$(k + 3)"] = attr(domain=cols[3], anchor="y$(k + 3)", range=[0, 500], tickvals=collect(TES_TICKS),
                                       ticktext="+" .* string.(TES_TICKS), showticklabels=bottom, ticks="")
        layout["yaxis$(k + 1)"] = attr(domain=y_dom, anchor="x$(k + 1)", showticklabels=true, ticks="",
                                       tickvals=tick_pos[r][1], ticktext=tick_pos[r][2],
                                       title=attr(text=var_label, standoff=4), automargin=true)
        layout["yaxis$(k + 2)"] = attr(domain=y_dom, anchor="x$(k + 2)", showticklabels=false, ticks="")
        layout["yaxis$(k + 3)"] = attr(domain=y_dom, anchor="x$(k + 3)", showticklabels=false, ticks="")
        layout["xaxis$(k + 4)"] = fig_axis(domain=cols[4], anchor="y$(k + 4)", showticklabels=bottom, tickangle=45,
                                           title=bottom ? var_label : "", showgrid=false)
        layout["yaxis$(k + 4)"] = fig_axis(domain=y_dom, anchor="x$(k + 4)", title=attr(text="z", standoff=4))
        append!(shapes, guide_lines("x$(k + 1)", "y$(k + 1)", TSS_TICKS))
        append!(shapes, guide_lines("x$(k + 2)", "y$(k + 2)", BODY_TICKS))
        append!(shapes, guide_lines("x$(k + 3)", "y$(k + 3)", TES_TICKS))
        push!(annotations, attr(text="<b>$('A' + r - 1)</b>   $(MARK_NAMES[m])", x=0, y=y_dom[2], xref="paper",
                                yref="paper", xanchor="left", yanchor="bottom", showarrow=false, xshift=-120))
    end
    layout[:shapes] = shapes
    layout[:annotations] = annotations
    return plot(data, layout)
end

# Human H3K9me3: heatmaps on top, box plot below
function human_figure(traces)
    tss, body, tes, boxes = traces
    cols = [[0.0, 0.22], [0.23, 0.77], [0.78, 1.0]]
    data = GenericTrace[]
    layout = fig_layout(margin=attr(l=130, r=20, t=20, b=90))
    row_ticks = nothing
    for (j, t) in enumerate((tss[1], body[1], tes[1]))
        t2 = copy_trace(t)
        row_ticks = upsample_rows!(t2)
        t2[:xaxis], t2[:yaxis] = "x$j", "y$j"
        t2[:colorscale] = HEATMAP_COLORS
        t2[:showscale] = j == 3
        t2[:colorbar] = attr(x=1.01, xanchor="left", y=0.78, yanchor="middle", len=0.45, thickness=12,
                             title=attr(text="z", side="right"), outlinewidth=0)
        push!(data, t2)
    end
    b = copy_trace(boxes[1])
    b[:xaxis], b[:yaxis] = "x4", "y4"
    b[:marker] = attr(size=3, color=PURPLE)
    b[:line] = attr(color=PURPLE, width=1.5)
    push!(data, b)
    y_heat = [0.55, 1.0]
    layout["xaxis1"] = attr(domain=cols[1], anchor="y1", range=[-500, 0], tickvals=collect(TSS_TICKS), ticktext=string.(TSS_TICKS), ticks="")
    layout["xaxis2"] = attr(domain=cols[2], anchor="y2", tickvals=collect(BODY_TICKS), ticktext=string.(BODY_TICKS) .* "%", ticks="")
    layout["xaxis3"] = attr(domain=cols[3], anchor="y3", range=[0, 500], tickvals=collect(TES_TICKS), ticktext="+" .* string.(TES_TICKS), ticks="")
    layout["yaxis1"] = attr(domain=y_heat, anchor="x1", ticks="", tickvals=row_ticks[1], ticktext=row_ticks[2],
                            title=attr(text="<i>dS</i>", standoff=4), automargin=true)
    layout["yaxis2"] = attr(domain=y_heat, anchor="x2", showticklabels=false)
    layout["yaxis3"] = attr(domain=y_heat, anchor="x3", showticklabels=false)
    layout["xaxis4"] = fig_axis(domain=[0, 1], anchor="y4", tickangle=45, title="<i>dS</i>", showgrid=false)
    layout["yaxis4"] = fig_axis(domain=[0, 0.38], anchor="x4", title="z")
    layout[:shapes] = vcat(guide_lines("x1", "y1", TSS_TICKS), guide_lines("x2", "y2", BODY_TICKS), guide_lines("x3", "y3", TES_TICKS))
    return plot(data, layout)
end

stars(p) = p < 0.001 ? "***" : p < 0.01 ? "**" : "*"

# Figure 3: pair-average expression by dS decile, with significant Dunn's comparisons as brackets
# (red: lower in the more diverged decile, green: higher)
function expression_ds_figure()
    pairs = CSV.read(joinpath(data_dir, "filtered/paralog_expr_ds.csv"), DataFrame)
    deciles = CSV.read(joinpath(STATS_DIR, "expression_vs_ds_deciles.csv"), DataFrame)
    dunn = CSV.read(joinpath(STATS_DIR, "expression_by_ds_dunn.csv"), DataFrame)
    decile = levelcode.(cut(pairs.dS, 10))
    labels = [@sprintf("%.2f–%.2f", lo, hi) for (lo, hi) in zip(deciles.dS_min, deciles.dS_max)]
    data = GenericTrace[fig_box(pairs.AvgExpr[decile .== d], labels[d], BLUE) for d in 1:10]

    sig = filter(row -> row.p_adj < 0.05, dunn)
    brackets = [(parse.(Int, split(c, " - "))..., p) for (c, p) in zip(sig.comparison, sig.p_adj)]
    sort!(brackets, by=b -> (b[2] - b[1], b[1]))
    levels = Vector{Vector{Tuple{Int, Int}}}()
    top = maximum(pairs.AvgExpr)
    step = 0.06 * (top - minimum(pairs.AvgExpr))
    shapes, annotations = [], []
    for (i, j, p) in brackets
        level = findfirst(l -> all(b -> j < b[1] || i > b[2], l), levels)
        isnothing(level) && (push!(levels, Tuple{Int, Int}[]); level = length(levels))
        push!(levels[level], (i, j))
        y = top + step * level
        color = deciles.MedianExpr[j] < deciles.MedianExpr[i] ? RED : GREEN
        push!(shapes, attr(type="path", path="M $(i - 1),$(y - step / 3) L $(i - 1),$y L $(j - 1),$y L $(j - 1),$(y - step / 3)",
                           xref="x", yref="y", line=attr(color=color, width=1.5)))
        push!(annotations, attr(text=stars(p), x=(i + j) / 2 - 1, y=y, xref="x", yref="y", yanchor="bottom",
                                showarrow=false, yshift=-4))
    end
    layout = fig_layout(xaxis=fig_axis(title="<i>dS</i> decile", tickangle=45, showgrid=false),
                        yaxis=fig_axis(title="log<sub>2</sub>(TPM + 0.5)",
                                       range=[minimum(pairs.AvgExpr) - step, top + step * (length(levels) + 1)]),
                        shapes=shapes, annotations=annotations)
    return plot(data, layout)
end

# Figure S13: dS distribution of the analysed pairs
function ds_histogram_figure()
    dicty = filter(row -> row.dS <= 3, CSV.read(joinpath(data_dir, "filtered/paralog_filt.tsv"), DataFrame)).dS
    human = filter(row -> row.dS <= 3, CSV.read(joinpath(data_dir, "filtered/human_paralog_info_filt.csv"), DataFrame)).dS
    data = [histogram(x=dicty, xbins=attr(start=0, size=0.1), marker_color=BLUE, xaxis="x1", yaxis="y1"),
            histogram(x=human, xbins=attr(start=0, size=0.02), marker_color=BLUE, xaxis="x2", yaxis="y2")]
    layout = fig_layout(xaxis1=fig_axis(anchor="y1", title="<i>dS</i>", domain=[0, 1]),
                        yaxis1=fig_axis(anchor="x1", domain=[0.58, 1], title="Pairs"),
                        xaxis2=fig_axis(anchor="y2", title="<i>dS</i>", domain=[0, 1]),
                        yaxis2=fig_axis(anchor="x2", domain=[0, 0.42], title="Pairs"),
                        bargap=0.05,
                        annotations=[attr(text="<i>D. discoideum</i>", x=0.98, y=1, xref="paper", yref="paper", xanchor="right", yanchor="top", showarrow=false),
                                     attr(text="<i>H. sapiens</i>", x=0.98, y=0.42, xref="paper", yref="paper", xanchor="right", yanchor="top", showarrow=false)])
    return plot(data, layout)
end

expr_label = "ln(TPM + 0.5)"
ds_label = "<i>dS</i>"

# Figures 1 and 2 (all stages) and S3-S5 / S7-S9 (fruiting body, mound, vegetative)
save_figure(enrichment_figure(load_traces("", "expr", "_All"), PANEL_ORDER, expr_label, PURPLE), "fig1"; height=1300)
save_figure(enrichment_figure(load_traces("", "dS", "_All"), PANEL_ORDER, ds_label, GREEN), "fig2"; height=1300)
for (name_expr, name_ds, stage) in [("figS03", "figS07", "F"), ("figS04", "figS08", "M"), ("figS05", "figS09", "V")]
    marks = stage == "V" ? PANEL_ORDER[1:3] : PANEL_ORDER  # no vegetative ATAC sample
    height = stage == "V" ? 1000 : 1300
    save_figure(enrichment_figure(load_traces("", "expr", "_$stage"), marks, expr_label, PURPLE), name_expr; height=height)
    save_figure(enrichment_figure(load_traces("", "dS", "_$stage"), marks, ds_label, GREEN), name_ds; height=height)
end

save_figure(expression_ds_figure(), "fig3"; width=800, height=650)
save_figure(human_figure(load_traces("human_", "dS", "")), "fig5"; width=800, height=750)
save_figure(ds_histogram_figure(), "figS13"; width=700, height=700)
