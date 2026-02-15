include("prelude.jl")

using Serialization
using PlotlyJS

"""Disable colorbars for all heatmap traces in a PlotlyJS figure.

Works with both `Plot` objects and displayed `SyncPlot` wrappers.
"""
function remove_heatmap_colorbars!(fig)
    # PlotlyJS sometimes wraps plots in `SyncPlot` for display; the underlying plot is in `.plot`.
    plt = hasproperty(fig, :plot) ? getproperty(fig, :plot) : fig

    # Hide per-trace colorbars.
    if hasproperty(plt, :data)
        for tr in plt.data
            ty = try
                tr[:type]
            catch
                nothing
            end
            if ty == "heatmap"
                tr[:showscale] = false
            end
        end
    end

    # Hide shared coloraxis colorbars too (coloraxis, coloraxis2, ...).
    if hasproperty(plt, :layout)
        if haskey(plt.layout, :coloraxis)
            plt.layout[:coloraxis][:showscale] = false
        end
        try
            for k in keys(plt.layout)
                if startswith(String(k), "coloraxis")
                    ca = plt.layout[k]
                    try
                        ca[:showscale] = false
                    catch
                    end
                end
            end
        catch
        end
    end

    return fig
end

save_plots = true
box_plots = true
tissue_type = "All" # Options: "All", "V", "M", "F"
x_tick_font = attr(family="Times New Roman", size=26)
y_tick_font = attr(family="Times New Roman", size=36)
title_font_size = 30
vertical_spacing = 0.1
x1_vals = [collect(-400:100:-100)...]
x2_vals = [collect(20:20:80)...]
x3_vals = [collect(100:100:400)...]
line_color = "lightgray"
line_type = "10px,30px"

ser_data_dir = "../../dicty_data/julia_serialized/"
plot_save_dir = "../../dicty_data/saved_figures/"

# Expression enrichment plots:
tss_expr_enrich = deserialize(joinpath(ser_data_dir, "tss_enrich_plots_expr_$tissue_type.jls"))
body_expr_enrich = deserialize(joinpath(ser_data_dir, "body_enrich_plots_expr_$tissue_type.jls"))
tes_expr_enrich_ends = deserialize(joinpath(ser_data_dir, "tes_enrich_plots_expr_$tissue_type.jls"))
bar_expr = deserialize(joinpath(ser_data_dir, "bar_plots_expr_$tissue_type.jls"))

# dS enrichment plots:
tss_dS_enrich = deserialize(joinpath(ser_data_dir, "tss_enrich_plots_dS_$tissue_type.jls"))
body_dS_enrich = deserialize(joinpath(ser_data_dir, "body_enrich_plots_dS_$tissue_type.jls"))
tes_dS_enrich_ends = deserialize(joinpath(ser_data_dir, "tes_enrich_plots_dS_$tissue_type.jls"))
bar_dS = deserialize(joinpath(ser_data_dir, "bar_plots_dS_$tissue_type.jls"))
tss_dS_enrich_human = deserialize(joinpath(ser_data_dir, "human_tss_enrich_plots_dS.jls"))
body_dS_enrich_human = deserialize(joinpath(ser_data_dir, "human_body_enrich_plots_dS.jls"))
tes_dS_enrich_human = deserialize(joinpath(ser_data_dir, "human_tes_enrich_plots_dS.jls"))
bar_dS_human = deserialize(joinpath(ser_data_dir, "human_bar_plots_dS.jls"))

sample_names = tissue_type == "V" ? ["H3K27ac", "H3K4me3", "H3K9me3"] : ["H3K4me3", "H3K27ac", "H3K9me3", "ATAC"]
for i in eachindex(sample_names)
    fig = make_subplots(
        rows=4, cols=3,
        column_widths=[0.5*0.75, 1*0.75, 0.5*0.75],
        # NOTE: Plotly applies row_heights bottom→top.
        # We want: (top) heatmap, bar, heatmap, bar (bottom)
        row_heights=[0.30, 0.30, 0.30, 0.30],
        horizontal_spacing=0.002,
        vertical_spacing=vertical_spacing,
        # shared_xaxes=true,
        # shared_yaxes=true,
        specs=[Spec(kind="heatmap") Spec(kind="heatmap") Spec(kind="heatmap");
               Spec(kind="bar", colspan=3) missing missing;
               Spec(kind="heatmap") Spec(kind="heatmap") Spec(kind="heatmap");
               Spec(kind="bar", colspan=3) missing missing]
    );

    # Expression enrichment plots:
    add_trace!(fig, tss_expr_enrich[1][i], row=1, col=1)
    add_trace!(fig, body_expr_enrich[1][i], row=1, col=2)
    add_trace!(fig, tes_expr_enrich_ends[1][i], row=1, col=3)
    add_trace!(fig, box_plots ? bar_expr[1][i] : bar_expr[i], row=2, col=1)
    relayout!(fig, yaxis=attr(tickfont=y_tick_font,
                              title="log(TPM + 0.5)",
                              titlefont=attr(family = "Times New Roman",
                              size=title_font_size)),
                xaxis=attr(range=[-500,0],
                            showticklabels=false),
                yaxis2=attr(showticklabels=false),
                xaxis2=attr(showticklabels=false),
                yaxis3=attr(showticklabels=false),
                xaxis3=attr(range=[0,500], 
                            showticklabels=false),
                yaxis4=attr(showticklabels=true,
                            tickfont=y_tick_font,
                            title="",
                            titlefont=attr(family="Times New Roman", size=title_font_size),
                            showline=true,
                            linecolor="black",
                            linewidth=2,
                            mirror=true),
                xaxis4=attr(tickfont=x_tick_font,
                            title="log(TPM + 0.5)",
                            titlefont=attr(family="Times New Roman", size=title_font_size),
                            showline=true,
                            linecolor="black",
                            linewidth=2,
                            mirror=true))

    # dS enrichment plots:
    add_trace!(fig, tss_dS_enrich[1][i], row=3, col=1)
    add_trace!(fig, body_dS_enrich[1][i], row=3, col=2)
    add_trace!(fig, tes_dS_enrich_ends[1][i], row=3, col=3)
    add_trace!(fig, box_plots ? bar_dS[1][i] : bar_dS[i], row=4, col=1)
    relayout!(fig, 
                yaxis5=attr(tickfont=y_tick_font, 
                               title="𝑑𝑆", 
                               titlefont=attr(family = "Times New Roman", 
                               size=title_font_size, 
                               weight="bold")),
                xaxis5=attr(range=[-500,0], 
                            tickmode="array", 
                            tickvals=x1_vals, 
                            ticktext=[string.(x1_vals)...], 
                            tickfont=x_tick_font),
                yaxis6=attr(showticklabels=false),
                xaxis6=attr(tickmode="array", 
                            tickvals=x2_vals, 
                            ticktext=string.(x2_vals) .* "%", 
                            tickfont=x_tick_font),
                yaxis7=attr(showticklabels=false),
                xaxis7=attr(range=[0,500], 
                            tickmode="array", 
                            tickvals=x3_vals, 
                            ticktext=["+" .* string.(x3_vals)...], 
                            tickfont=x_tick_font),
                yaxis8=attr(showticklabels=true,
                            tickfont=y_tick_font,
                            title="",
                            titlefont=attr(family="Times New Roman", size=title_font_size),
                            showline=true,
                            linecolor="black",
                            linewidth=2,
                            mirror=true),
                xaxis8=attr(tickfont=x_tick_font,
                            title="𝑑𝑆",
                            titlefont=attr(family="Times New Roman", size=title_font_size, weight="bold"),
                            showline=true,
                            linecolor="black",
                            linewidth=2,
                            mirror=true),
                title=sample_names[i] * " $tissue_type", 
                titlefont=attr(size=title_font_size, 
                               family="Times New Roman"))

    relayout!(fig, showlegend=false,
                   plot_bgcolor="rgba(0,0,0,0)",
                   height=2000)
    
        relayout!(fig, shapes = reduce(vcat, 
        [
        [   # Vertical lines for TSS heatmaps (rows 1 & 3, col 1)
        (type = "line",
         x0 = val, 
         x1 = val, 
         y0 = 0, 
         y1 = 1, 
         xref = "x", 
         yref = "y domain", 
         line = (color = line_color, 
             dash = line_type))
        for val in x1_vals],
        [
        (type = "line",
         x0 = val, 
         x1 = val, 
         y0 = 0, 
         y1 = 1, 
         xref = "x5", 
         yref = "y5 domain", 
         line = (color = line_color, 
             dash = line_type))
        for val in x1_vals],
        [   # Vertical lines for body heatmaps (rows 1 & 3, col 2)
        (type = "line",
         x0 = val, 
         x1 = val, 
         y0 = 0, 
         y1 = 1, 
         xref = "x2", 
         yref = "y2 domain", 
         line = (color = line_color, 
             dash = line_type))
        for val in x2_vals],
        [
        (type = "line",
         x0 = val, 
         x1 = val, 
         y0 = 0, 
         y1 = 1, 
         xref = "x6", 
         yref = "y6 domain", 
         line = (color = line_color, 
             dash = line_type))
        for val in x2_vals],
        [   # Vertical lines for TES heatmaps (rows 1 & 3, col 3)
        (type = "line",
         x0 = val,
         x1 = val,
         y0 = 0,
         y1 = 1,
         xref = "x3",
         yref = "y3 domain",
         line = (color = line_color, 
             dash = line_type))
        for val in x3_vals],
        [
        (type = "line",
         x0 = val,
         x1 = val,
         y0 = 0,
         y1 = 1,
         xref = "x7",
         yref = "y7 domain",
         line = (color = line_color, 
             dash = line_type))
        for val in x3_vals]
        ]
        ))

    remove_heatmap_colorbars!(fig)
    
    if save_plots
        savefig(fig, joinpath(plot_save_dir, "combined_plot_$(sample_names[i])_dS_$tissue_type.html"))
    
    end
    display(fig)
end

#________________________________________________________________________________________
# Human:
fig_h = make_subplots(
    rows=2, cols=3,
    column_widths=[0.5*0.75, 1*0.75, 0.5*0.75],
    # NOTE: Plotly applies row_heights bottom→top.
    row_heights=[0.30, 0.30],
    horizontal_spacing=0.005,
    vertical_spacing=vertical_spacing,
        # shared_xaxes=true,
        # shared_yaxes=true,
    specs=[Spec(kind="heatmap") Spec(kind="heatmap") Spec(kind="heatmap");
           Spec(kind="bar", colspan=3) missing missing]
);

# dS enrichment plots:
add_trace!(fig_h, tss_dS_enrich_human[1][1], row=1, col=1)
add_trace!(fig_h, body_dS_enrich_human[1][1], row=1, col=2)
add_trace!(fig_h, tes_dS_enrich_human[1][1], row=1, col=3)
add_trace!(fig_h, box_plots ? bar_dS_human[1][1] : bar_dS_human[1], row=2, col=1)
relayout!(fig_h, yaxis=attr(tickfont=y_tick_font, 
                            title="𝑑𝑆", 
                            titlefont=attr(family = "Times New Roman", 
                            size=title_font_size)),
            xaxis=attr(range=[-500,0],
                        tickmode="array", 
                        tickvals=x1_vals, 
                        ticktext=[string.(x1_vals)...], 
                        tickfont=x_tick_font),
            yaxis2=attr(showticklabels=false),
            xaxis2=attr(tickmode="array", 
                        tickvals=x2_vals, 
                        ticktext=string.(x2_vals) .* "%", 
                        tickfont=x_tick_font),
            yaxis3=attr(showticklabels=false),
            xaxis3=attr(range=[0,500], 
                        tickmode="array", 
                        tickvals=x3_vals, 
                        ticktext=["+" .* string.(x3_vals)...], 
                        tickfont=x_tick_font),
            yaxis4=attr(showticklabels=true,
                        tickfont=y_tick_font,
                        title="",
                        titlefont=attr(family="Times New Roman", size=title_font_size),
                        showline=true,
                        linecolor="black",
                        linewidth=2,
                        mirror=true),
            xaxis4=attr(tickfont=x_tick_font,
                        title="𝑑𝑆",
                        titlefont=attr(family="Times New Roman", size=title_font_size, weight="bold"),
                        showline=true,
                        linecolor="black",
                        linewidth=2,
                        mirror=true),
            title="H3K9me3",
            titlefont=attr(size=title_font_size, 
                           family="Times New Roman"))

relayout!(fig_h, showlegend=false)
relayout!(fig_h, title="H3K9me3")
relayout!(fig_h, plot_bgcolor="rgba(0,0,0,0)", height=1200)
relayout!(fig_h, shapes = reduce(vcat, 
    [
        [   # Vertical lines for subplot 1
        (type = "line",
         x0 = val, 
         x1 = val, 
         y0 = 0, 
         y1 = 1, 
         xref = "x", 
         yref = "y domain", 
         line = (color = line_color, 
                 dash = line_type))
        for val in x1_vals],
        [   # Vertical lines for subplot 2
        (type = "line",
         x0 = val,
         x1 = val,
         y0 = 0,
         y1 = 1,
         xref = "x2",
         yref = "y2 domain",
         line = (color = line_color,
                 dash = line_type))
        for val in x2_vals],
        [   # Vertical lines for subplot 3
        (type = "line",
         x0 = val, 
         x1 = val, 
         y0 = 0, 
         y1 = 1, 
         xref = "x3", 
         yref = "y3 domain", 
         line = (color = line_color,
                 dash = line_type))
        for val in x3_vals]
    ]
))

remove_heatmap_colorbars!(fig_h)

if save_plots
    savefig(fig_h, joinpath(plot_save_dir, "combined_plot_H3K9me3_dS_human.html"))

end
display(fig_h)