using CSV
using DataFrames
using StatsBase
using PlotlyJS
using RollingFunctions
using DSP

function cor_range(ds_vals, id_vals, ds_threshold, cor_func::F=cor) where F <: Function
    supr = findfirst(ds_vals .> ds_threshold) - 1
    last_ind = isnothing(supr) ? length(ds_vals) : supr - 1
    
    if last_ind == 0
        error("no dS vals less than or equal to threshold")
    
    end
    
    ds_vals_slice = ds_vals[1:last_ind]
    id_vals_slice = id_vals[1:last_ind]

    return cor_func(ds_vals_slice, id_vals_slice)
end

# Paralog file
full_paralog_file = "../../../../data/AX4/genome_ver_2_7/biomart/genes_v52/d_disc_paralogs_biomart_ensembl_protist_ver_52.txt"
filt_paralog_file = "./data/"

# Load paralog data
paralog_data = CSV.read(full_paralog_file, DataFrame, delim='\t')
ds_estimates = CSV.read("./data/filtered/dS_df.csv", DataFrame)

select!(paralog_data, [1, 15, 4, 5, 6, 7])
rename!(paralog_data, ["GeneID", "ParalogID", "dS","IDqt","IDtq","Type"])

# Remove gene splits and non-paralogs
filter!(row -> !any(ismissing.(Vector(row))) && row.Type != "gene_split", paralog_data)
# filter!(row -> mean([row.IDqt, row.IDtq]) >= 30, paralog_data)

# Sort by dS and calculate correlations
sort!(paralog_data, :dS)
ds_vals = Vector{Float64}(paralog_data.dS)
id_vals = map(pair -> max(pair...), zip(paralog_data.IDqt,paralog_data.IDtq))
roll_cor = rolling(cor, ds_vals, id_vals, 1000)
id_windowed = rollmean(id_vals, 1000)
ds_windowed = rollmean(ds_vals, 1000)

mean_filter = DSP.gaussian(10000,1)
mean_filter = mean_filter ./ sum(mean_filter)

plot(scatter(x=filtfilt(mean_filter, filtfilt(mean_filter, ds_vals))[1:5:end], y=filtfilt(mean_filter, filtfilt(mean_filter, id_vals))[1:5:end], mode="markers"), Layout(xaxis=attr(range=[0,8])))
plot(scatter(x=filtfilt(mean_filter, id_windowed), y=filtfilt(mean_filter, roll_cor), mode="markers"))
plot(scatter(x=filtfilt(mean_filter, ds_windowed), y=filtfilt(mean_filter, roll_cor), mode="markers"))