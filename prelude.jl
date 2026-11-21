# Set working directory
using Pkg
cd(@__DIR__)
Pkg.activate(".")

# Custom library
using BioinfoTools
using .LoadGFF
using .GenomicData
using .GenomeTypes

# Common packages
using CSV
using DataFrames
using HypothesisTests
using MultipleTesting
using PlotlyJS
using Random
using StatsBase
using Statistics

# Constants
z_min = -1
z_max = 1

# Fixed seed so permutation-test p-values are reproducible
Random.seed!(1)

# Command-line options shared by the analysis scripts:
#   first positional argument: life-cycle stage ("All", "V", "M" or "F")
#   --yn00: use YN00 dS (old_pairs_with_new_ds_dicty.csv) instead of Ensembl Compara dS
stage_arg(default="All") = (pos = filter(a -> !startswith(a, "--"), ARGS); isempty(pos) ? default : pos[1])
use_yn00() = "--yn00" in ARGS

# Statistical results are written here as CSV files
const STATS_DIR = mkpath("../../dicty_data/stats")
write_stats(name, df) = CSV.write(joinpath(STATS_DIR, "$name.csv"), df)

# Effect size r = |Z| / √N for a Mann–Whitney U test (normal approximation; for the exact test, σ is
# computed without a tie correction).
function mwu_r(test)
    nx, ny = test.nx, test.ny
    σ = hasproperty(test, :sigma) ? test.sigma : sqrt(nx * ny * (nx + ny + 1) / 12)
    return abs(test.U - nx * ny / 2) / σ / √(nx + ny)
end

# One row summarising a Mann–Whitney U test of x vs y
mwu_row(label, x, y) = (test = MannWhitneyUTest(x, y);
    (Test = label, n1 = length(x), n2 = length(y), median1 = median(x), median2 = median(y),
     U = test.U, r = mwu_r(test), p = pvalue(test)))

# One row summarising a Kruskal–Wallis test, with η² = (H - k + 1) / (N - k)
function kw_row(label, test::KruskalWallisTest)
    k, N = length(test.n_i), sum(test.n_i)
    return (Test = label, H = test.H, df = test.df, N = N, eta2 = (test.H - k + 1) / (N - k), p = pvalue(test))
end

# Figure style and export
include("figure_style.jl")
