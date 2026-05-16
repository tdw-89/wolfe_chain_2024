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
using StatsBase
using Statistics

# Constants
z_min = -1
z_max = 1