using Serialization
using PlotlyJS

save_plots = true
box_plots = true
x_tick_font = attr(family="Times New Roman", size=26)
y_tick_font = attr(family="Times New Roman", size=26)
title_font_size = 30
x1_vals = [collect(-400:100:-100)...]
x2_vals = [collect(20:20:80)...]
x3_vals = [collect(100:100:400)...]
line_color = "lightgray"
line_type = "10px,30px"

ser_data_dir = "./data/julia_serialized/"
plot_save_dir = "./data/saved_figures/"

# Expression enrichment plots:
tss_expr_enrich = deserialize(joinpath(ser_data_dir, "tss_enrich_plots_expr.jls"))
body_expr_enrich = deserialize(joinpath(ser_data_dir, "body_enrich_plots_expr.jls"))
tes_expr_enrich_ends = deserialize(joinpath(ser_data_dir, "tes_enrich_plots_expr.jls"))
bar_expr = deserialize(joinpath(ser_data_dir, "bar_plots_expr.jls"))

# dS enrichment plots:
tss_dS_enrich = deserialize(joinpath(ser_data_dir, "tss_enrich_plots_dS.jls"))
body_dS_enrich = deserialize(joinpath(ser_data_dir, "body_enrich_plots_dS.jls"))
tes_dS_enrich_ends = deserialize(joinpath(ser_data_dir, "tes_enrich_plots_dS.jls"))
bar_dS = deserialize(joinpath(ser_data_dir, "bar_plots_dS.jls"))
tss_dS_enrich_human = deserialize(joinpath(ser_data_dir, "human_tss_enrich_plots_dS.jls"))
body_dS_enrich_human = deserialize(joinpath(ser_data_dir, "human_body_enrich_plots_dS.jls"))
tes_dS_enrich_human = deserialize(joinpath(ser_data_dir, "human_tes_enrich_plots_dS.jls"))
bar_dS_human = deserialize(joinpath(ser_data_dir, "human_bar_plots_dS.jls"))

sample_names = ["K27ac", 
                "K4me3", 
                "K9me3", 
                "ATAC"]
for i in eachindex(sample_names)
    fig = make_subplots(
        rows=2, cols=4,
        column_widths=[0.5*0.75, 1*0.75, 0.5*0.75, 0.25],
        row_heights=[0.5, 0.5],
        horizontal_spacing=0.002,
        vertical_spacing=0.02,
        # shared_xaxes=true,
        # shared_yaxes=true,
        specs=[Spec(kind="heatmap") Spec(kind="heatmap") Spec(kind="heatmap") Spec(kind="bar");
               Spec(kind="heatmap") Spec(kind="heatmap") Spec(kind="heatmap") Spec(kind="bar")]
    );

    # Expression enrichment plots:
    add_trace!(fig, tss_expr_enrich[1][i], row=1, col=1)
    add_trace!(fig, body_expr_enrich[1][i], row=1, col=2)
    add_trace!(fig, tes_expr_enrich_ends[1][i], row=1, col=3)
    add_trace!(fig, box_plots ? bar_expr[i].plot.data[1] : bar_expr[i], row=1, col=4)
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
                yaxis4=attr(showticklabels=false),
                xaxis4=attr(range=[0,4],
                            showticklabels=false))

    # dS enrichment plots:
    add_trace!(fig, tss_dS_enrich[1][i], row=2, col=1)
    add_trace!(fig, body_dS_enrich[1][i], row=2, col=2)
    add_trace!(fig, tes_dS_enrich_ends[1][i], row=2, col=3)
    add_trace!(fig, box_plots ? bar_dS[i].plot.data[1] : bar_dS[i], row=2, col=4)
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
                yaxis8=attr(showticklabels=false),
                xaxis8=attr(range=[0,4], tickfont=x_tick_font),
                title=sample_names[i], 
                titlefont=attr(size=title_font_size, 
                               family="Times New Roman"))

    relayout!(fig, showlegend=false,
                   plot_bgcolor="rgba(0,0,0,0)")
    
    relayout!(fig, shapes = reduce(vcat, 
    [
        [   # Vertical lines for subplots 1 & 4
        (type = "line",
         x0 = val, 
         x1 = val, 
         y0 = 0, 
         y1 = 1, 
         xref = "x", 
         yref = "paper", 
         line = (color = line_color, 
                 dash = line_type))
        for val in x1_vals],
        [   # Vertical lines for subplots 2 & 5
        (type = "line",
         x0 = val, 
         x1 = val, 
         y0 = 0, 
         y1 = 1, 
         xref = "x2", 
         yref = "paper", 
         line = (color = line_color, 
                 dash = line_type))
        for val in x2_vals],
        [   # Vertical lines for subplots 3 & 6
        (type = "line",
         x0 = val,
         x1 = val,
         y0 = 0,
         y1 = 1,
         xref = "x3",
         yref = "paper",
         line = (color = line_color, 
                 dash = line_type))
        for val in x3_vals]
    ]
    ))
    
    if save_plots
        savefig(fig, joinpath(plot_save_dir, "proposal_combined_plot_$(sample_names[i])_dS.html"))
    
    end
    display(fig)
end

# Human:
fig_h = make_subplots(
        rows=1, cols=4,
        column_widths=[0.5*0.75, 1*0.75, 0.5*0.75, 0.25],
        horizontal_spacing=0.005,
        # shared_xaxes=true,
        # shared_yaxes=true,
        specs=[Spec(kind="heatmap") Spec(kind="heatmap") Spec(kind="heatmap") Spec(kind="bar")]
);

# dS enrichment plots:
add_trace!(fig_h, tss_dS_enrich_human[1][1], row=1, col=1)
add_trace!(fig_h, body_dS_enrich_human[1][1], row=1, col=2)
add_trace!(fig_h, tes_dS_enrich_human[1][1], row=1, col=3)
add_trace!(fig_h, box_plots ? bar_dS_human[1].plot.data[1] : bar_dS_human[1], row=1, col=4)
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
            yaxis4=attr(showticklabels=false),
            xaxis4=attr(range=[0,4], tickfont=x_tick_font),
            title="H3K9me3",
            titlefont=attr(size=title_font_size, 
                           family="Times New Roman"))

relayout!(fig_h, showlegend=false)
relayout!(fig_h, title="H3K9me3")
relayout!(fig_h, plot_bgcolor="rgba(0,0,0,0)")
relayout!(fig_h, shapes = reduce(vcat, 
    [
        [   # Vertical lines for subplots 1 & 4
        (type = "line",
         x0 = val, 
         x1 = val, 
         y0 = 0, 
         y1 = 1, 
         xref = "x", 
         yref = "paper", 
         line = (color = line_color, 
                 dash = line_type))
        for val in x1_vals],
        [   # Vertical lines for subplots 2 & 5
        (type = "line",
         x0 = val,
         x1 = val,
         y0 = 0,
         y1 = 1,
         xref = "x2",
         yref = "paper",
         line = (color = line_color,
                 dash = line_type))
        for val in x2_vals],
        [   # Vertical lines for subplots 3 & 6
        (type = "line",
         x0 = val, 
         x1 = val, 
         y0 = 0, 
         y1 = 1, 
         xref = "x3", 
         yref = "paper", 
         line = (color = line_color,
                 dash = line_type))
        for val in x3_vals]
    ]
))

if save_plots
    savefig(fig_h, joinpath(plot_save_dir, "combined_plot_H3K9me3_dS_human.html"))

end
display(fig_h)
