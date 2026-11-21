# Shared style and export for the manuscript figures (included by prelude.jl).
#
# Figures are written as vector PDF to ../../dicty_data/figures/<name>.pdf for \includegraphics. Sizes are in
# CSS pixels; a figure FIG_WIDTH px wide placed at \textwidth (~6.5 in) gets ~9 pt text from FIG_FONT_SIZE.

const FIG_DIR = mkpath("../../dicty_data/figures")
const FIG_WIDTH = 1000
const FIG_FONT_SIZE = 18
const FIG_FONT = attr(family="Times New Roman", size=FIG_FONT_SIZE, color="black")
# Plasma (plotly.js has no built-in name for it)
const HEATMAP_COLORS = [[0.0, "#0d0887"], [0.111, "#46039f"], [0.222, "#7201a8"], [0.333, "#9c179e"], [0.444, "#bd3786"],
                        [0.556, "#d8576b"], [0.667, "#ed7953"], [0.778, "#fb9f3a"], [0.889, "#fdca26"], [1.0, "#f0f921"]]
const BLUE, RED, GREEN, PURPLE = "#3b5bdb", "#e03131", "#2f9e44", "#9c6ade"

const KALEIDO_RESTARTED = Ref(false)

function save_figure(fig, name; width=FIG_WIDTH, height=600)
    # PlotlyJS starts Kaleido with a bundled MathJax that leaves a "failed to load" banner in exported
    # figures; restart it once without MathJax
    if !KALEIDO_RESTARTED[]
        PlotlyJS.PlotlyKaleido.restart(mathjax=false)
        KALEIDO_RESTARTED[] = true
    end
    savefig(fig, joinpath(FIG_DIR, "$name.pdf"); width=width, height=height)
end

# Axis with a black frame and light grid
fig_axis(; kwargs...) = attr(; merge((showline=true, linecolor="black", linewidth=1, mirror=true, ticks="outside",
                                      gridcolor="#e6e6e6", zeroline=false, automargin=true), kwargs)...)

fig_layout(; kwargs...) = Layout(; merge((font=FIG_FONT, plot_bgcolor="white", paper_bgcolor="white", showlegend=false,
                                          margin=attr(l=70, r=20, t=20, b=70)), kwargs)...)

fig_box(y, name, color; kwargs...) = box(; y=y, name=name, marker_color=color, line=attr(color=color, width=1.5),
                                         boxpoints="outliers", marker_size=3, kwargs...)
