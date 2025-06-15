#= 
    This module contains a set of functions useful for calculating and plotting
    enrichment of peaks (from ChIP, ATAC, etc.) or calculating peak coverage
=#

using Printf
using Pipe
using PlotlyJS
using CategoricalArrays
using NaturalSort
using Interpolations
using HypothesisTests
using Random
using Statistics

if !@isdefined ChromData

    include("genomic_data.jl")
end

if !@isdefined FLOAT_RE

    const FLOAT_RE = r"[0-9]+\.[0-9]+"
    const LBOUND_RE = r"[\[(]{1}"
    const UBOUND_RE = r"[\])]{1}"
    const QUANT_RANGE_RE = r"^.*(?<lbound>[\[(]{1})\s*(?<lval>[\-0-9]+\.[0-9]+)\s*,\s*(?<uval>[\-0-9]+\.[0-9]+)\s*(?<ubound>[\])]{1})$"
end

const RANGE_PRECISION = 2  # Precision for rounding quantile range values

include("misc_utils.jl")

function parse_quantile(q::String; digs::Int=RANGE_PRECISION)
    mat = match(QUANT_RANGE_RE, q)
    if mat === nothing
        @warn "Warning: label '$q' does not match the expected quantile range format. Keeping it as is."
        return q
    end

    lbound = mat["lbound"]
    lval = round(parse(Float64, mat["lval"]), digits=digs)
    uval = round(parse(Float64, mat["uval"]), digits=digs)
    ubound = mat["ubound"]
    return lbound * string(lval) * ", " * string(uval) * ubound

end

parse_quantile(quant_labels::Vector{String}; digs::Int=3) = [parse_quantile(q, digs=digs) for q in quant_labels]

function sortNparsequantrange(quant_labels::Vector{String}; digs::Int=RANGE_PRECISION)

    sorted_labels = sort(quant_labels, lt = natural)
    ret_labels = String[]

    for label in sorted_labels

        mat = match(QUANT_RANGE_RE, label)
        if mat === nothing
            push!(ret_labels, label)  # If it doesn't match the range pattern, show a warning
            @warn "Warning: label '$label' does not match the expected quantile range format. Keeping it as is."
            continue
        end

        lbound = mat["lbound"]
        ubound = mat["ubound"]
        lval = round(parse(Float64, mat["lval"]), digits=digs)
        uval = round(parse(Float64, mat["uval"]), digits=digs)
        temp_str = string(lbound, lval, ", ", uval, ubound)

        push!(ret_labels, temp_str)
    end

    return ret_labels
end

"""
NOTES: 
    Currently this method treats the first exon as the region TSS:TSS+500
    This method will throw an error if:
        * the first exon is less than 500 bp when selecting 'first_exon' as the region,
        * the upstream region is <= 250 bp when selecting 'upstream' as the region,
        * the upstream region is <= 250 bp or the gene body is less than 150 bp if 'promoter' as the region,
    
"""
function getrange(gene::Gene, range_type::String, use_binsignals::Bool=true)
    
    temp_gene_start = gene.gene_start
    temp_gene_end = gene.gene_end
    temp_region_start = gene.regions[1].region_start

    if use_binsignals

        if isempty(gene.regions[1].binsignals)

            error("missing 'binsignal' data for gene $(gene.id)")
        else

            sig_length = length(gene.regions[1].binsignals[1])
        end
    else

        if isempty(gene.regions[1].signals)

            error("missing 'signal' data for gene $(gene.id)")
        else

            sig_length = length(gene.regions[1].signals[1])
        end
    end

    if range_type == "upstream"

        if gene.strand == '+'

            if (temp_gene_start - temp_region_start) - 250 > 0

                return 1:((temp_gene_start - temp_region_start) - 250)
            else

                error("gene $(gene.id) has an upstream region less than 250 bp")
            end
        else

            if (temp_gene_end - temp_region_start) + 250 <= sig_length

                return ((temp_gene_end - temp_region_start) + 250):sig_length
            else
                    
                error("gene $(gene.id) has an upstream region less than 250 bp")
            end
        end

    elseif range_type == "downstream"

        if gene.strand == '+'

            return (temp_gene_end - temp_region_start):sig_length
        else

            return 1:(temp_gene_start - temp_region_start)
        end

    elseif range_type == "gene_body"

        
        return (temp_gene_start - temp_region_start):(temp_gene_end - temp_region_start)
        

    elseif range_type == "promoter"

        if gene.strand == '+'

            reg_start = (temp_gene_start - temp_region_start) - 250
            reg_end = (temp_gene_start - temp_region_start) + 150

            if reg_start > 0 && reg_end <= sig_length

                return reg_start:reg_end
            else

                if reg_start <= 0

                    error("gene $(gene.id) has a promoter region less than or equal to 250 bp")
                else
                        
                    error("gene $(gene.id) has a gene body less than 150 bp")
                end
            end

        else

            reg_end = (temp_gene_end - temp_region_start) + 250
            reg_start = (temp_gene_end - temp_region_start) - 150

            if reg_end <= sig_length && reg_start > 0

                return reg_start:reg_end
            else

                if reg_end > sig_length

                    error("gene $(gene.id) has a promoter region less than 250 bp")
                else
                        
                    error("gene $(gene.id) has a gene body less than or equal to 150 bp")
                end
            end
        end 

    elseif range_type == "first_exon"

        if gene.strand == '+'

            reg_end = (temp_gene_start - temp_region_start) + 500
            reg_start = (temp_gene_start - temp_region_start)
            
            if reg_end <= sig_length

                return reg_start:reg_end
            else

                error("gene $(gene.id) has a first exon less than 500 bp")
            end
        else

            reg_end = (temp_gene_end - temp_region_start)
            reg_start = (temp_gene_end - temp_region_start) - 500

            if reg_start > 0

                return reg_start:reg_end
            else

                error("gene $(gene.id) has a first exon less than 500 bp")
            end
        end
    else
        error("Invalid peak range string. Must be one of \"downstream\", \"upstream\", \"promoter\", \"gene_body\", or \"first_exon\".")
    end
end

"""
Get the signal for a gene in a given sample in the given region around the TSS specified by 'upstream_lim' and 'downstream_lim'.
    The upstream and downstream limits are automatically adjusted to account for strandedness, so provide those arguments as if the gene is on the '+' strand.
    The default values for 'upstream_lim' and 'downstream_lim' are 1999 and 2000, respectively, which will return a 4000 bp region centered on the TSS.
        The return value is a vector with the signal values for the given region presented as if in the + strand orientation.
            e.g. 5' -> 3' => 1 -> 4000
"""
function getsiginrange_old(gene::Gene, gene_range::GeneRange, sample_ind=1; signal_type::String="peak")

    if sample_ind ∉ 1:length(gene.samples)

        error("Invalid sample index for gene $(gene.id).")
    end

    # Get the signal vector
    if signal_type == "peak"

        if isempty(gene.regions[1].binsignals)

            error("No peak data for gene $(gene.id) in sample $sample_ind.")
        end

        signal = gene.regions[1].binsignals[sample_ind]
    else

        error("not implemented")
    end

    if gene.strand == '+'
        
        if typeof(gene_range.range_start) == REGION

            range_start = 1 + gene_range.start_offset
        elseif typeof(gene_range.range_start) == TSS

            range_start = gene.gene_start - gene.regions[1].region_start + 1 + gene_range.start_offset
        elseif typeof(gene_range.range_start) == TES

            range_start = gene.gene_end - gene.regions[1].region_start + 1 + gene_range.start_offset
        else

            error("unrecognized range start $(gene_range.range_start)")
        end

        if typeof(gene_range.range_stop) == REGION

            range_stop = gene.regions[1].region_end - gene.regions[1].region_start + 1 + gene_range.stop_offset
        elseif typeof(gene_range.range_stop) == TSS

            range_stop = gene.gene_start - gene.regions[1].region_start + 1 + gene_range.stop_offset
        elseif typeof(gene_range.range_stop) == TES
                
            range_stop = gene.gene_end - gene.regions[1].region_start + 1 + gene_range.stop_offset
        else

            error("unrecognized range stop $(gene_range.range_stop)")
        end

        try
            return signal[range_start:range_stop]
        catch BoundsError
            println("region $gene_range for gene $(gene.id) out of bounds.")
            return missing
        end
    elseif gene.strand == '-'

        if typeof(gene_range.range_start) == REGION
            range_start = gene.regions[1].region_end - gene.regions[1].region_start + 1 - gene_range.start_offset

        elseif typeof(gene_range.range_start) == TSS
            range_start = gene.gene_end - gene.regions[1].region_start + 1 - gene_range.start_offset

        elseif typeof(gene_range.range_start) == TES
            range_start = gene.gene_start - gene.regions[1].region_start + 1 - gene_range.start_offset

        else

            error("unrecognized range start $(gene_range.range_start)")
        end

        if typeof(gene_range.range_stop) == REGION
            range_stop = 1 - gene_range.stop_offset

        elseif typeof(gene_range.range_stop) == TSS
            range_stop = gene.gene_end - gene.regions[1].region_start + 1 - gene_range.stop_offset

        elseif typeof(gene_range.range_stop) == TES
            range_stop = gene.gene_start - gene.regions[1].region_start + 1 - gene_range.stop_offset

        else

            error("unrecognized range stop $(gene_range.range_stop)")
        end

        try
            return reverse(signal[range_stop:range_start])
        catch BoundsError
            println("region $gene_range for gene $(gene.id) out of bounds.")
            return missing
        end
    else

        error("invalid strand $(gene.strand)")
    end
end

function siginrange_old(gene::Gene, sig_range::GeneRange, sample_ind=1; peak_data::Bool=true)
    if sample_ind ∉ 1:length(gene.samples)

        error("Invalid sample index for gene $(gene.id).")
    end

    # Get the signal vector
    if peak_data == true

        if isempty(gene.regions[1].binsignals)

            error("No peak data for gene $(gene.id) in sample $sample_ind.")
        end

        signal = gene.regions[1].binsignals[sample_ind]
    else

        error("not implemented")
    end

    if gene.strand == '+'
        
        if typeof(sig_range.range_start) == REGION
            range_start = 1 + sig_range.start_offset

        elseif isa(sig_range.range_start, TSS) && isa(sig_range.start_offset_type, PERCENTAGE)
            gene_len = gene.gene_end - gene.gene_start + 1
            gene_len_frac = min(round(Int, gene_len / 100 * sig_range.start_offset), gene_len)
            range_start = gene.gene_start - gene.regions[1].region_start + 1 + gene_len_frac

        elseif isa(sig_range.range_start, TSS)
            range_start = gene.gene_start - gene.regions[1].region_start + 1 + sig_range.start_offset

        elseif isa(sig_range.range_start, TES)
            range_start = gene.gene_end - gene.regions[1].region_start + 1 + sig_range.start_offset

        else

            error("unrecognized range start $(sig_range.range_start)")
        end

        if isa(sig_range.range_stop, REGION)

            range_stop = gene.regions[1].region_end - gene.regions[1].region_start + 1 + sig_range.stop_offset
        elseif isa(sig_range.range_stop, TSS) && isa(sig_range.stop_offset_type, PERCENTAGE)
            gene_len = gene.gene_end - gene.gene_start + 1
            gene_len_frac = min(round(Int, gene_len / 100 * sig_range.start_offset), gene_len)
            range_stop = gene.gene_start - gene.regions[1].region_start + 1 + gene_len_frac
        elseif isa(sig_range.range_stop, TSS)

            range_stop = gene.gene_start - gene.regions[1].region_start + 1 + sig_range.stop_offset
        elseif isa(sig_range.range_stop, TES)
                
            range_stop = gene.gene_end - gene.regions[1].region_start + 1 + sig_range.stop_offset
        else

            error("unrecognized range stop $(sig_range.range_stop)")
        end

        try 
            sig = signal[range_start:range_stop]
            return true
        catch BoundsError
            return false
        end
    elseif gene.strand == '-'

        if isa(sig_range.range_start, REGION)
            range_start = gene.regions[1].region_end - gene.regions[1].region_start + 1 - sig_range.start_offset

        elseif isa(sig_range.range_start, TSS) && isa(sig_range.stop_offset_type, PERCENTAGE)

            range_start = gene.gene_end - gene.regions[1].region_start + 1 - sig_range.start_offset

        elseif isa(sig_range.range_start, TSS)
            range_start = gene.gene_end - gene.regions[1].region_start + 1 - sig_range.start_offset

        elseif isa(sig_range.range_start, TES)
            range_start = gene.gene_start - gene.regions[1].region_start + 1 - sig_range.start_offset

        else
            error("unrecognized range start $(sig_range.range_start)")

        end

        if typeof(sig_range.range_stop) == REGION
            range_stop = 1 - sig_range.stop_offset

        elseif typeof(sig_range.range_stop) == TSS
            range_stop = gene.gene_end - gene.regions[1].region_start + 1 - sig_range.stop_offset

        elseif typeof(sig_range.range_stop) == TES
            range_stop = gene.gene_start - gene.regions[1].region_start + 1 - sig_range.stop_offset

        else

            error("unrecognized range stop $(sig_range.range_stop)")
        end

        try
            reverse(signal[range_stop:range_start])
            return true
        catch err
            if isa(err, BoundsError)
                return false
            else
                throw(err)
            end
        end
    else

        error("invalid strand $(gene.strand)")
    end
end

function getsiginrange(gene::Gene, sig_range::GeneRange, sample_ind=1; peak_data::Bool=true, clamped::Bool=true, pad_clamped::Bool=true)
    if sample_ind ∉ 1:length(gene.samples)

        error("Invalid sample index for gene $(gene.id).")
    end

    # Get the signal vector
    if !peak_data
        error("not implemented")

    end

    if isempty(gene.regions[1].binsignals)
        error("No peak data for gene $(gene.id) in sample $sample_ind.")

    end

    signal = gene.strand == '+' ? gene.regions[1].binsignals[sample_ind] : reverse(gene.regions[1].binsignals[sample_ind])
    gene_len = gene.gene_end - gene.gene_start + 1
    tss = gene.strand == '+' ? (gene.gene_start - gene.regions[1].region_start) + 1 : (gene.regions[1].region_end - gene.gene_end) + 1
    tes = tss + gene_len - 1
    range_start = isa(sig_range.range_start, REGION) ? 1 : 
                   isa(sig_range.range_start, TES) ? tes + sig_range.start_offset :
                   isa(sig_range.range_start, TSS) && isa(sig_range.start_offset_type, PERCENTAGE) ? tss + round(Int, (gene_len / 100) * sig_range.start_offset) :
                   tss + sig_range.start_offset
    range_stop = isa(sig_range.range_stop, REGION) ? length(signal) : 
                 isa(sig_range.range_stop, TES) ? tes + sig_range.stop_offset :
                 isa(sig_range.range_stop, TSS) && isa(sig_range.stop_offset_type, PERCENTAGE) ? tss + round(Int, (gene_len / 100) * sig_range.stop_offset) :
                 tss + sig_range.stop_offset

    pad_left = 0
    pad_right = 0
    if clamped
        if pad_clamped
            start_clamped = clamp(range_start, 1, length(signal) + 1)
            stop_clamped = clamp(range_stop, 0, length(signal))
            # println("$tss $tes $range_start $range_stop $start_clamped $stop_clamped")
            pad_left = abs(start_clamped - range_start)
            pad_right = abs(stop_clamped - range_stop)
            range_start = start_clamped
            range_stop = stop_clamped
        else
            range_start = clamp(range_start, 1, length(signal) + 1)
            range_stop = clamp(range_stop, 0, length(signal))
        end
    end
    try
        return vcat(
            ones(pad_left) * -1, signal[range_start:range_stop], ones(pad_right) * -1)
    catch err
        if isa(err, BoundsError)
            return missing
        else
            throw(err)
        end
    end
end

function siginrange(gene::Gene, sig_range::GeneRange, sample_ind=1; peak_data::Bool=true)
    if sample_ind ∉ 1:length(gene.samples)

        error("Invalid sample index for gene $(gene.id).")
    end

    # Get the signal vector
    if !peak_data
        error("not implemented")

    end

    if isempty(gene.regions[1].binsignals)
        error("No peak data for gene $(gene.id) in sample $sample_ind.")

    end

    signal = gene.strand == '+' ? gene.regions[1].binsignals[sample_ind] : reverse(gene.regions[1].binsignals[sample_ind])
    gene_len = gene.gene_end - gene.gene_start + 1
    tss = gene.strand == '+' ? (gene.gene_start - gene.regions[1].region_start) + 1 : (gene.regions[1].region_end - gene.gene_end) + 1
    tes = gene.strand == '+' ? (gene.gene_end - gene.regions[1].region_start) + 1 : (gene.regions[1].region_end - gene.gene_start) + 1
    range_start = isa(sig_range.range_start, REGION) ? 1 + sig_range.start_offset : 
                   isa(sig_range.range_start, TES) ? tes + sig_range.start_offset :
                   isa(sig_range.range_start, TSS) && isa(sig_range.start_offset_type, PERCENTAGE) ? tss + round(Int, gene_len / sig_range.start_offset * 100) :
                   tss + sig_range.start_offset
    range_stop = isa(sig_range.range_stop, REGION) ? length(signal) + sig_range.stop_offset : 
                 isa(sig_range.range_stop, TES) ? tes + sig_range.stop_offset :
                 isa(sig_range.range_stop, TSS) && isa(sig_range.stop_offset_type, PERCENTAGE) ? tss + round(Int, gene_len / sig_range.stop_offset * 100) :
                 tss + sig_range.stop_offset

    try 
        sig = signal[range_start:range_stop]
        return true
    catch err
        if isa(err, BoundsError)
            return false
        else
            throw(err)
        end
    end
end

function to_percent(signal_vec::Vector{R}) where R <: Union{Real}

    signal_interp = interpolate((1:length(signal_vec),), signal_vec, Gridded(Linear()))
    new_inds = range(1, length(signal_vec), 100)
    return signal_interp[new_inds]
end

"""
    Average the peak signal for each gene across each sample (in the specified gene region).
    The indices for the different samples for each peak type need to be specified, along
    with the peak type name.
    NOTES: 
        * Sample index tuples must correspond 1:1 to the sample type names
        * This function only calculates combinatorial presence in a specific range of each gene
        * This function filters out genes that have no peaks of any kind in any sample
"""
function findcombinations(qualified_genes::Vector{Gene}, sample_inds::Vector{T}, peak_type_names::Vector{String}, peak_range::Union{R, String}; print_only::Bool=false) where {T <: Tuple{Vararg{Int}}, R <: AbstractRange} 
   
    n_peak_types = length(sample_inds)
    n_genes = length(qualified_genes)
    peak_name_combos = collect(powerset(peak_type_names, 1))
    peak_name_combos = map(x -> reduce(*, x), peak_name_combos)
    peak_num_combos = collect(powerset(1:n_peak_types, 1))
    overlap_df = hcat(DataFrame("GeneID" => [gene.id for gene in qualified_genes]), DataFrame([combo => BitVector(zeros(n_genes)) for combo in peak_name_combos]))
    temp_enrich_vec = BitVector(zeros(n_peak_types))
    to_exclude = Int[]

    if typeof(peak_range) == String
    
        lowercase(peak_range) in ["downstream", "upstream", "promoter", "gene_body", "first_exon"] ? nothing : error("Invalid peak range string. Must be one of \"downstream\", \"upstream\", \"promoter\", \"gene_body\", or \"first_exon\".") 
        
        for i in eachindex(qualified_genes)
        
            temp_range = getrange(qualified_genes[i], peak_range)

            for n in 1:n_peak_types

                for k in sample_inds[n]
    
                    temp_enrich_vec[n] = sum(qualified_genes[i].regions[1].binsignals[k][temp_range]) > 0                    
                end

                # temp_enrich_vec[n] /= length(sample_inds[n])
            end
            
            for j in 1:(size(overlap_df, 2) - 1)

                overlap_df[i, j+1] = reduce(&, temp_enrich_vec[peak_num_combos[j]])
            end

            temp_enrich_vec .= 0
        end
    else

        for i in eachindex(qualified_genes)
        
            for n in 1:n_peak_types

                for k in sample_inds[n]

                    if qualified_genes[i].strand == '+'
                    
                        temp_enrich_vec[n] = sum(qualified_genes[i].regions[1].binsignals[k][peak_range]) > 0 
                    else

                        temp_enrich_vec[n] = sum(reverse(qualified_genes[i].regions[1].binsignals[k])[peak_range]) > 0 
                    end
                end

                # temp_enrich_vec[n] /= length(sample_inds[n])
            end
            
            for j in 1:(size(overlap_df, 2) - 1)

                overlap_df[i, j+1] = reduce(&, temp_enrich_vec[peak_num_combos[j]])
            end

            temp_enrich_vec .= 0
        end
    end

    if print_only

        comb_perc_med = map(x -> median(x[x .!= 0]), eachcol(overlap_df)[2:end])
        comb_perc_mean = map(x -> mean(x[x .!= 0]), eachcol(overlap_df)[2:end])

        for i in 1:length(comb_perc_med)
        
            println("$(names(overlap_df)[i+1]) median enrichment: $(comb_perc_med[i])")
            println("$(names(overlap_df)[i+1]) mean enrichment: $(comb_perc_mean[i])")
        end
    else

        return overlap_df[Not(to_exclude),:]
    end
end

function combination_prop(gene_list::Vector{Gene}, 
    comb_ind_tuples::Vector{T}, 
    region::Union{String, R}, 
    sample_names::Vector{String}) where {T <: Tuple, R <: UnitRange}    

    comb_df = findcombinations(gene_list, comb_ind_tuples, sample_names, region)
    n_genes = length(gene_list)
    combination_prop = map(col -> sum(col)/length(col), eachcol(comb_df)[2:end])
    combination_names = names(comb_df)[2:end]
    comb_df_summary_all = DataFrame("CombName" => [combination_names[i] for i in eachindex(combination_names)], "CombProp" => [combination_prop[i] for i in eachindex(combination_prop)])
    comb_df_summary_all = vcat(comb_df_summary_all, DataFrame("CombName" => "n genes", "CombProp" => n_genes))
    return comb_df_summary_all    
end

function plot_enrich_region(paralog_df::DataFrame, gene_list::Vector{Gene}, sample_groups::Vector{T}, group_regions::Vector{GeneRange}; fold_change_over_mean::Bool=false, 
                                                                        global_means::Union{Float64, Vector{Float64}, Nothing}=nothing,
                                                                        z_min::Int=0, z_max::Int=4,
                                                                        return_figs::Bool=false,
                                                                        save_plots::Bool=false,
                                                                        ind_var_col=3) where T <: Union{Tuple, Vector{Int}}

    fig_vec = [GenericTrace[], Layout[]]
    for (i, sample_inds) in enumerate(sample_groups) 

        n_quantiles = 10
        
        sample_name = gene_list[1].samples[sample_inds[1]]

        gene_ds_quantile_labels = String[]
        gene_ds_quantiles = Int[]
        quantile_legend = String[]

        while isempty(gene_ds_quantile_labels) && n_quantiles > 1

            try

                gene_ds_quantile_labels = cut(paralog_df[!,ind_var_col], n_quantiles)
                gene_ds_quantiles = levelcode.(gene_ds_quantile_labels)
                quantile_legend = sortNparsequantrange(unique(gene_ds_quantile_labels))
                
                
            catch ArgumentError

                n_quantiles -= 1
                gene_ds_quantile_labels = String[]
                gene_ds_quantiles = Int[]
                quantile_legend = String[]
            end
        end

        if isempty(gene_ds_quantile_labels)

            println("ERROR: No quantiles could be created for $sample_name.")
            continue
        end

        n_quantiles = length(unique(gene_ds_quantiles))
        x_range = to_vector(group_regions[i])
        positional_count_mat = zeros(length(x_range), n_quantiles)

        for d in 1:n_quantiles

            # Find the indices of all the gene pairs in the paralog df whose ds is in quantile 'd'.
            # Stop after identifying ⌊n pairs / 2⌋ pairs and their associated indices in the 
            # 'gene_list'. 
            gene_pair_inds = findall(gene_ds_quantiles .== d)
            # shuffle!(gene_pair_inds)
            
            qual_gene_ind_pairs = Vector{Tuple{Int, Int}}()
            # max_pairs = length(gene_pair_inds)
            for ind in gene_pair_inds

                gene_1 = paralog_df[ind, 1]
                gene_2 = paralog_df[ind, 2]

                qual_gene_1_ind = findfirst([gene.id == gene_1 for gene in gene_list])
                qual_gene_2_ind = findfirst([gene.id == gene_2 for gene in gene_list])

                if !(isnothing(qual_gene_1_ind)) && !(isnothing(qual_gene_2_ind))

                    push!(qual_gene_ind_pairs, (qual_gene_1_ind, qual_gene_2_ind))
                end

                # if length(qual_gene_ind_pairs) >= max_pairs
                    
                #     break
                # end
            end

            for g in eachindex(qual_gene_ind_pairs)
                temp_gene = gene_list[qual_gene_ind_pairs[g][1]]
                temp_paralog = gene_list[qual_gene_ind_pairs[g][2]]
                gene_pos_count = [getsiginrange(temp_gene, group_regions[i], sample_ind) for sample_ind in sample_inds]
                paralog_pos_count = [getsiginrange(temp_paralog, group_regions[i], sample_ind) for sample_ind in sample_inds]
                if any(ismissing.(gene_pos_count)) || any(ismissing.(paralog_pos_count))
                    @warn "skipping gene pair $(temp_gene.id) - $(temp_paralog.id) because of missing data in specified range."
                    continue
                end

                gene_pos_count = mean(reduce(hcat, gene_pos_count), dims=2)[:,1]
                paralog_pos_count = mean(reduce(hcat, paralog_pos_count), dims=2)[:,1]
                positional_count_mat[:,d] .= positional_count_mat[:,d] .+ mean([gene_pos_count paralog_pos_count], dims=2)[:,1]

            end

            positional_count_mat[:,d] .= positional_count_mat[:,d] ./ length(qual_gene_ind_pairs)
        end

        if fold_change_over_mean
                
            if isnothing(global_means)
                
                mean_val = mean(positional_count_mat)
            elseif typeof(global_means) == Vector{Float64}
                
                mean_val = global_means[i]
            else

                mean_val = global_means
            end

            positional_count_mat = positional_count_mat / mean_val
        end

        

        layout_1 = Layout(
            xaxis_side="bottom",
            title="$sample_name positional enrichment vs. ds $n_quantiles-quantile",
            width=1500,
            height=750,

        )

        if save_plots

            savefig(plot(heatmap(x=x_range, y=quantile_legend, z=positional_count_mat', zmin=z_min, zmax=z_max), layout_1), joinpath(plot_save_dir, "$(sample_name)_positional_enrichment_vs_ds_quantile_$(n_quantiles)_quantiles_heatmap.html"))
        end

        if return_figs
            push!(fig_vec[1], 
                heatmap(
                    x=x_range,
                    y=quantile_legend,
                    z=positional_count_mat',
                    zmin=z_min,
                    zmax=z_max,
                    # colorscale=[[0, "rgb(255,255,255)"], [1, "rgb(255,0,0)"]]),
                ))
            push!(fig_vec[2], layout_1)

        else
            display(
                plot(
                    heatmap(
                        x=x_range,
                        y=quantile_legend,
                        z=positional_count_mat',
                        zmin=z_min,
                        zmax=z_max,
                        # colorscale=[[0, "rgb(255,255,255)"], [1, "rgb(255,0,0)"]]),
                    ), layout_1
                )
            )

        end
    end

    if return_figs
        return fig_vec

    end
end

function plot_enrich_percent(paralog_df::DataFrame, gene_list::Vector{Gene}, sample_groups::Vector{T}; fold_change_over_mean::Bool=false, 
                                                                        global_means::Union{Float64, Vector{Float64}, Nothing}=nothing,
                                                                        z_min::Int=0, z_max::Int=4,
                                                                        return_figs::Bool=false,
                                                                        save_plots::Bool=false,
                                                                        ind_var_col=3) where T <: Union{Tuple, Vector{Int}}

    fig_vec = [GenericTrace[], Layout[]]
    for (i, sample_inds) in enumerate(sample_groups) 

        n_quantiles = 10
        
        sample_name = gene_list[1].samples[sample_inds[1]]

        gene_ds_quantile_labels = String[]
        gene_ds_quantiles = Int[]
        quantile_legend = String[]

        while isempty(gene_ds_quantile_labels) && n_quantiles > 1

            try

                gene_ds_quantile_labels = cut(paralog_df[!,ind_var_col], n_quantiles)
                gene_ds_quantiles = levelcode.(gene_ds_quantile_labels)
                quantile_legend = sortNparsequantrange(unique(gene_ds_quantile_labels))
                
            catch ArgumentError

                n_quantiles -= 1
                gene_ds_quantile_labels = String[]
                gene_ds_quantiles = Int[]
                quantile_legend = String[]
            end
        end

        if isempty(gene_ds_quantile_labels)

            println("ERROR: No quantiles could be created for $sample_name.")
            continue
        end

        n_quantiles = length(unique(gene_ds_quantiles))
        x_range = 1:100
        positional_count_mat = zeros(length(x_range), n_quantiles)

        for d in 1:n_quantiles

            # Find the indices of all the gene pairs in the paralog df whose ds is in quantile 'd'.
            # Stop after identifying ⌊n pairs / 2⌋ pairs and their associated indices in the 
            # 'gene_list'. 
            gene_pair_inds = findall(gene_ds_quantiles .== d)
            # shuffle!(gene_pair_inds)
            
            qual_gene_ind_pairs = Vector{Tuple{Int, Int}}()
            # max_pairs = length(gene_pair_inds)
            for ind in gene_pair_inds

                gene_1 = paralog_df[ind, 1]
                gene_2 = paralog_df[ind, 2]

                qual_gene_1_ind = findfirst([gene.id == gene_1 for gene in gene_list])
                qual_gene_2_ind = findfirst([gene.id == gene_2 for gene in gene_list])

                if !(isnothing(qual_gene_1_ind)) && !(isnothing(qual_gene_2_ind))

                    push!(qual_gene_ind_pairs, (qual_gene_1_ind, qual_gene_2_ind))
                end

                # if length(qual_gene_ind_pairs) >= max_pairs
                    
                #     break
                # end
            end

            for g in eachindex(qual_gene_ind_pairs)
                temp_gene = gene_list[qual_gene_ind_pairs[g][1]]
                temp_paralog = gene_list[qual_gene_ind_pairs[g][2]]
                gene_pos_count = [getsiginrange(temp_gene, GeneRange(TSS(), TES(), 0, 0), sample_ind) for sample_ind in sample_inds]
                paralog_pos_count = [getsiginrange(temp_paralog, GeneRange(TSS(), TES(), 0, 0), sample_ind) for sample_ind in sample_inds]
                if any(ismissing.(gene_pos_count)) || any(ismissing.(paralog_pos_count))
                    @warn "skipping gene pair $(temp_gene.id) - $(temp_paralog.id) because of missing data in specified range."
                    continue
                end

                gene_pos_count = to_percent(mean(reduce(hcat, gene_pos_count), dims=2)[:,1])
                paralog_pos_count = to_percent(mean(reduce(hcat, paralog_pos_count), dims=2)[:,1])
                positional_count_mat[:,d] .= positional_count_mat[:,d] .+ mean([gene_pos_count paralog_pos_count], dims=2)[:,1]

            end

            positional_count_mat[:,d] .= positional_count_mat[:,d] ./ length(qual_gene_ind_pairs)
        end

        if fold_change_over_mean
                
            if isnothing(global_means)
                
                mean_val = mean(positional_count_mat)
            elseif typeof(global_means) == Vector{Float64}
                
                mean_val = global_means[i]
            else

                mean_val = global_means
            end

            positional_count_mat = positional_count_mat / mean_val
        end

        

        layout_1 = Layout(
            xaxis_side="bottom",
            title="$sample_name positional enrichment vs. ds $n_quantiles-quantile",
            width=1500,
            height=750,

        )

        if save_plots

            savefig(plot(heatmap(x=x_range, y=quantile_legend, z=positional_count_mat', zmin=z_min, zmax=z_max), layout_1), joinpath(plot_save_dir, "$(sample_name)_positional_enrichment_vs_ds_quantile_$(n_quantiles)_quantiles_heatmap.html"))
        end

        if return_figs
            push!(fig_vec[1], 
                heatmap(
                    x=x_range,
                    y=quantile_legend,
                    z=positional_count_mat',
                    zmin=z_min,
                    zmax=z_max,
                    # colorscale=[[0, "rgb(255,255,255)"], [1, "rgb(255,0,0)"]]),
                ))
            push!(fig_vec[2], layout_1)

        else
            display(
                plot(
                    heatmap(
                        x=x_range,
                        y=quantile_legend,
                        z=positional_count_mat',
                        zmin=z_min,
                        zmax=z_max,
                        # colorscale=[[0, "rgb(255,255,255)"], [1, "rgb(255,0,0)"]]),
                    ), layout_1
                )
            )

        end
    end

    if return_figs
        return fig_vec

    end
end

function plot_enrich_expr_region(expr_df::DataFrame, gene_list::Vector{Gene}, sample_groups::Vector{T}, group_regions::Vector{GeneRange}; fold_change_over_mean::Bool=false, 
                                                                        global_means::Union{Float64, Vector{Float64}, Nothing}=nothing,
                                                                        z_min::Int=0, z_max::Int=4,
                                                                        return_figs::Bool=false,
                                                                        save_plots::Bool=false) where T <: Union{Tuple, Vector{Int}}

    fig_vec = [GenericTrace[], Layout[]]
    for (i, sample_inds) in enumerate(sample_groups) 

        sample_name = gene_list[1].samples[sample_inds[1]]
        n_quantiles = 10
        gene_expr_quantile_labels = String[]
        gene_expr_quantiles = Int[]
        quantile_legend = String[]

        while isempty(gene_expr_quantile_labels) && n_quantiles > 1

            try
                gene_expr_quantile_labels = cut(expr_df.Avg, n_quantiles)
                gene_expr_quantiles = levelcode.(gene_expr_quantile_labels)
                quantile_legend = sortNparsequantrange(unique(gene_expr_quantile_labels))
                
            catch ArgumentError
                n_quantiles -= 1
                gene_expr_quantile_labels = String[]
                gene_expr_quantiles = Int[]
                quantile_legend = String[]

            end
        end

        if isempty(gene_expr_quantile_labels)
            println("ERROR: No quantiles could be created for $sample_name.")
            continue

        end

        expr_df_quantiles = @pipe copy(expr_df) |> insertcols!(_, :quantile => gene_expr_quantiles)
        n_quantiles = length(unique(gene_expr_quantiles))
        
        x_range = to_vector(group_regions[i])
        positional_count_mat = zeros(length(x_range), n_quantiles)

        for d in 1:n_quantiles
            gene_inds = findall(expr_df_quantiles.quantile .== d)
            qual_genes = [gene for gene in gene_list if gene.id in expr_df_quantiles[gene_inds, :GeneID]]

            for gene in qual_genes
                gene_sample_counts = [getsiginrange(gene, group_regions[i], sample_ind) for sample_ind in sample_inds]
                if any(ismissing.(gene_sample_counts))
                    @warn "skipping gene $(gene.id) because it is missing data in one or more samples in the specified range."
                    continue
                end

                pair_pos_count = mean(reduce(hcat, gene_sample_counts), dims=2)[:,1]
                positional_count_mat[:,d] .= positional_count_mat[:,d] .+ pair_pos_count

            end
            positional_count_mat[:,d] .= positional_count_mat[:,d] ./ length(qual_genes)

        end

        if fold_change_over_mean
                
            if isnothing(global_means)
                mean_val = mean(positional_count_mat)

            elseif typeof(global_means) == Vector{Float64}
                mean_val = global_means[i]

            else
                mean_val = global_means

            end
            positional_count_mat = positional_count_mat / mean_val

        end

        layout_1 = Layout(
            xaxis_side="bottom",
            title="$sample_name positional enrichment vs. expr $n_quantiles-quantile",
            width=1500,
            height=750,

        )

        if save_plots
            savefig(plot(heatmap(x=x_range, y=quantile_legend, z=positional_count_mat', zmin=z_min, zmax=z_max), layout_1), joinpath(plot_save_dir, "$(sample_name)_positional_enrichment_vs_expr_quantile_$(n_quantiles)_quantiles_heatmap.html"))
        
        end

        if return_figs
            push!(fig_vec[1], 
                heatmap(
                    x=x_range,
                    y=quantile_legend,
                    z=positional_count_mat',
                    zmin=z_min,
                    zmax=z_max,
                    # colorscale=[[0, "rgb(255,255,255)"], [1, "rgb(255,0,0)"]]),
                ))
            push!(fig_vec[2], layout_1)

        else
            display(
                plot(
                    heatmap(
                        x=x_range,
                        y=quantile_legend,
                        z=positional_count_mat',
                        zmin=z_min,
                        zmax=z_max,
                        # colorscale=[[0, "rgb(255,255,255)"], [1, "rgb(255,0,0)"]]),
                    ), layout_1
                )
            )

        end
    end

    if return_figs
        return fig_vec

    end
end

function plot_enrich_expr_percent(expr_df::DataFrame, gene_list::Vector{Gene}, sample_groups::Vector{T}; fold_change_over_mean::Bool=false, 
                                    global_means::Union{Float64, Vector{Float64}, Nothing}=nothing,
                                    z_min::Int=0, z_max::Int=4,
                                    return_figs::Bool=false,
                                    save_plots::Bool=false) where T <: Union{Tuple, Vector{Int}}
    
    fig_vec = [GenericTrace[], Layout[]]
    for (i, sample_inds) in enumerate(sample_groups) 

        sample_name = gene_list[1].samples[sample_inds[1]]
        n_quantiles = 10
        gene_expr_quantile_labels = String[]
        gene_expr_quantiles = Int[]
        quantile_legend = String[]

        while isempty(gene_expr_quantile_labels) && n_quantiles > 1

            try
                gene_expr_quantile_labels = cut(expr_df.Avg, n_quantiles)
                gene_expr_quantiles = levelcode.(gene_expr_quantile_labels)
                quantile_legend = sortNparsequantrange(unique(gene_expr_quantile_labels))
                
            catch ArgumentError
                n_quantiles -= 1
                gene_expr_quantile_labels = String[]
                gene_expr_quantiles = Int[]
                quantile_legend = String[]

            end
        end

        if isempty(gene_expr_quantile_labels)
            println("ERROR: No quantiles could be created for $sample_name.")
            continue

        end

        expr_df_quantiles = @pipe copy(expr_df) |> insertcols!(_, :quantile => gene_expr_quantiles)
        n_quantiles = length(unique(gene_expr_quantiles))
        
        x_range = 1:100
        positional_count_mat = zeros(length(x_range), n_quantiles)

        for d in 1:n_quantiles
            gene_inds = findall(expr_df_quantiles.quantile .== d)
            qual_genes = [gene for gene in gene_list if gene.id in expr_df_quantiles[gene_inds, :GeneID]]

            for gene in qual_genes
                gene_sample_counts = [getsiginrange(gene, GeneRange(TSS(), TES(), 0, 0), sample_ind) for sample_ind in sample_inds]
                if any(ismissing.(gene_sample_counts))
                    @warn "skipping gene $(gene.id) because it is missing data in one or more samples in the specified range."
                    continue
                end

                pos_count = to_percent(mean(reduce(hcat, gene_sample_counts), dims=2)[:,1])
                if any(isnan.(pos_count)) # DEBUG REMOVE
                    return gene, sample_inds, gene_sample_counts # DEBUG REMOVE
                end # DEBUG REMOVE
                positional_count_mat[:,d] .= positional_count_mat[:,d] .+ pos_count

            end
            positional_count_mat[:,d] .= positional_count_mat[:,d] ./ length(qual_genes)

        end

        if fold_change_over_mean
                
            if isnothing(global_means)
                mean_val = mean(positional_count_mat)

            elseif typeof(global_means) == Vector{Float64}
                mean_val = global_means[i]

            else
                mean_val = global_means

            end
            positional_count_mat = positional_count_mat / mean_val

        end

        layout_1 = Layout(
            xaxis_side="bottom",
            title="$sample_name positional enrichment vs. expr $n_quantiles-quantile",
            width=1500,
            height=750,

        )

        if save_plots
            savefig(plot(heatmap(x=x_range, y=quantile_legend, z=positional_count_mat', zmin=z_min, zmax=z_max), layout_1), joinpath(plot_save_dir, "$(sample_name)_positional_enrichment_vs_expr_quantile_$(n_quantiles)_quantiles_heatmap.html"))
        
        end

        if return_figs
            push!(fig_vec[1], 
                heatmap(
                    x=x_range,
                    y=quantile_legend,
                    z=positional_count_mat',
                    zmin=z_min,
                    zmax=z_max,
                    # colorscale=[[0, "rgb(255,255,255)"], [1, "rgb(255,0,0)"]]),
                ))
            push!(fig_vec[2], layout_1)

        else
            display(
                plot(
                    heatmap(
                        x=x_range,
                        y=quantile_legend,
                        z=positional_count_mat',
                        zmin=z_min,
                        zmax=z_max,
                        # colorscale=[[0, "rgb(255,255,255)"], [1, "rgb(255,0,0)"]]),
                    ), layout_1
                )
            )

        end
    end

    if return_figs
        return fig_vec

    end
end

"""
    Plot the average signal for each gene pair in the specified region across the samples in `sample_inds`.
    The average signal is calculated as the log2 of the ratio of the average signal for each gene pair to the average signal given.
    The function will also plot the average signal for each quantile of ds.
    The function will return a vector of figures if return_figs is true, otherwise it will display the figures.
"""
function plot_bar(paralog_df::DataFrame, 
                  gene_list::Vector{Gene},
                  sample_ind_groups::Vector{T},
                  region_list::Vector{GeneRange},
                  global_means::Vector{Float64},
                  y_range::Union{Nothing, Vector{Int}},
                  return_figs::Bool=true,
                  horizontal::Bool=false,
                  save_plots::Bool=false,
                  box_plots::Bool=true,
                  ind_var_col=3) where T <: Union{Tuple, Vector{Int}}
    
    @assert length(region_list) == length(global_means) == length(sample_ind_groups) "length of 'region_list', 'global_means', and 'sample_ind_groups' must be equal."
    
    fig_vec = []

    p_values = []

    means_vecs = Vector{Vector{Vector{Float64}}}() # DEBUG REMOVE

    for (r, sample_inds) in enumerate(sample_ind_groups)

        n_quantiles = 10
        
        sample_name = gene_list[1].samples[sample_inds[1]]

        gene_ds_quantile_labels = String[]
        gene_ds_quantiles = Int[]
        quantile_legend = String[]

        while isempty(gene_ds_quantile_labels) && n_quantiles > 1

            try

                gene_ds_quantile_labels = cut(paralog_df[!,ind_var_col], n_quantiles)
                gene_ds_quantiles = levelcode.(gene_ds_quantile_labels)
                quantile_legend = sortNparsequantrange(unique(gene_ds_quantile_labels))
                
            catch ArgumentError

                n_quantiles -= 1
                gene_ds_quantile_labels = String[]
                gene_ds_quantiles = Int[]
                quantile_legend = String[]
            end
        end

        if isempty(gene_ds_quantile_labels)

            println("ERROR: No quantiles could be created for $sample_name.")
            continue
        end

        n_quantiles = length(unique(gene_ds_quantiles))
        
        pair_means = [Vector{Float64}() for i in 1:n_quantiles]

        for d in 1:n_quantiles

            # Find the indices of all the gene pairs in the paralog df whose ds is in quantile 'd'.
            # Stop after identifying ⌊n pairs / 2⌋ pairs and their associated indices in the 
            # 'gene_list'. 
            gene_pair_inds = findall(gene_ds_quantiles .== d)
            # shuffle!(gene_pair_inds)

            # q_1 = Tuple{String, String}[] # DEBUG REMOVE
            # q_2 = Tuple{String, String}[] # DEBUG REMOVE
            
            qual_gene_ind_pairs = Vector{Tuple{Int, Int}}()
            # max_pairs = length(gene_pair_inds)
            for ind in gene_pair_inds

                gene_1 = paralog_df[ind, 1]
                gene_2 = paralog_df[ind, 2]

                qual_gene_1_ind = findfirst([gene.id == gene_1 for gene in gene_list])
                qual_gene_2_ind = findfirst([gene.id == gene_2 for gene in gene_list])

                if !(isnothing(qual_gene_1_ind)) && !(isnothing(qual_gene_2_ind))

                    push!(qual_gene_ind_pairs, (qual_gene_1_ind, qual_gene_2_ind))
                end
            end

            for i in eachindex(qual_gene_ind_pairs)

                # Get the average count within the specified region for each sample in sample_inds. Then average across samples.
                gene_sigs = [getsiginrange(gene_list[qual_gene_ind_pairs[i][1]], region_list[r], sample_ind) for sample_ind in sample_inds]
                if !any(ismissing.(gene_sigs))
                    pair_pos_count_gene = mean(reduce(vcat, gene_sigs))

                end

                # Repeat the for the second gene in the paralog pair
                paralog_sigs = [getsiginrange(gene_list[qual_gene_ind_pairs[i][2]], region_list[r], sample_ind) for sample_ind in sample_inds]
                if !any(ismissing.(paralog_sigs))
                    pair_pos_count_paralog = mean(reduce(vcat, paralog_sigs))

                end

                # Add the average of the two to the vector of average counts for the current quantile
                if !any(ismissing.(gene_sigs)) && !any(ismissing.(paralog_sigs))
                    pair_average = mean([pair_pos_count_gene, pair_pos_count_paralog])
                    push!(pair_means[d], pair_average/global_means[r])
                
                else
                    @warn "skipping $(gene_list[qual_gene_ind_pairs[i][1]].id) and $(gene_list[qual_gene_ind_pairs[i][2]].id)
                    One or both genes do not have data in the specified region for samples $(gene_list[qual_gene_ind_pairs[i][2]].samples[sample_inds])."
                
                end
            end
        end

        quant_ranges = String[]

        for d in 1:n_quantiles
            gene_pair_inds = findall(gene_ds_quantiles .== d)
            quant_min = minimum(paralog_df[!,ind_var_col][gene_pair_inds])
            quant_max = maximum(paralog_df[!,ind_var_col][gene_pair_inds])
            
            if d == n_quantiles
            
                quant_range = @sprintf("[%.2f - %.2f]", quant_min, quant_max)
            else
            
                quant_range = @sprintf("[%.2f - %.2f)", quant_min, quant_max)
            end
            
            push!(quant_ranges, quant_range)
        end

        if box_plots
            vals = []
            quantile_nums = []
            for (i, mean_dist) in enumerate(pair_means)
                if !isempty(mean_dist)
                    vals = vcat(vals, mean_dist)
                    quantile_nums = vcat(quantile_nums, fill(i, length(mean_dist)))
                end
            end
            # Create a box plot for the values across quantiles
            temp_df = DataFrame(
                :X => quantile_nums,
                :Y => vals
            )
            if horizontal
                bar_plt = plot(temp_df, x=:Y, y=:X, kind="box", orientation = "h")
            else
                bar_plt = plot(temp_df, x=:X, y=:Y, kind="box")
            end

        elseif isnothing(y_range)
            bar_plt = bar(
                x=[mean(pair_mean_dist) for pair_mean_dist in pair_means],
                y=quant_ranges,
                orientation = horizontal ? "h" : "v"
            )
        elseif length(y_range) == 2
            if horizontal
                bar_plt = bar(
                    x=[mean(pair_mean_dist) for pair_mean_dist in pair_means],
                    y=quant_ranges,
                    orientation = "h",
                    xaxis=attr(range=y_range)
                )
            else
                bar_plt = bar(
                    x=[mean(pair_mean_dist) for pair_mean_dist in pair_means],
                    y=quant_ranges,
                    orientation = "v",
                    yaxis=attr(range=y_range)
                )
            end
        else
            @warn "Invalid y_range. Must be either nothing or a vector of length 2."
            bar_plt = bar(
                x=[mean(pair_mean_dist) for pair_mean_dist in pair_means],
                y=quant_ranges,
                orientation = horizontal ? "h" : "v"
            )
        end

        if save_plots
            savefig(box_plots ? bar_plt : plot(bar_plt), joinpath(plot_save_dir, "$(sample_name)_region_enrichment_vs_ds_$(n_quantiles)-quantile_average_bar.html"))

        end

        if return_figs
            push!(fig_vec, bar_plt)

        else
            display(plot(bar_plt))
        
        end

        push!(p_values, KruskalWallisTest(pair_means...))
        push!(means_vecs, pair_means)
    end

    if return_figs
        return fig_vec, p_values, means_vecs

    end

    return p_values
end

function plot_bar_expr(expr_df::DataFrame, 
                  gene_list::Vector{Gene}, 
                  sample_ind_groups::Vector{T}, 
                  region_list::Vector{GeneRange}, 
                  global_means::Vector{Float64}, 
                  y_range::Union{Nothing, Vector{Int}},
                  return_figs::Bool=false,
                  horizontal::Bool=false,
                  save_plots::Bool=false,
                  box_plots::Bool=true) where T <: Union{Tuple, Vector{Int}}
    
    @assert length(region_list) == length(global_means) == length(sample_ind_groups) "length of 'region_list', 'global_means', and 'sample_ind_groups' must be equal."

    if return_figs
        fig_vec = []

    end

    p_values = []

    means_vecs = Vector{Vector{Vector{Float64}}}()

    n_quantiles = 10
    gene_expr_quantile_labels = String[]
    gene_expr_quantiles = Int[]
    quantile_legend = String[]

    while isempty(gene_expr_quantile_labels) && n_quantiles > 1
        try
            gene_expr_quantile_labels = cut(expr_df.Avg, n_quantiles)
            gene_expr_quantiles = levelcode.(gene_expr_quantile_labels)
            quantile_legend = sortNparsequantrange(unique(gene_expr_quantile_labels))
        catch ArgumentError
            n_quantiles -= 1
            gene_expr_quantile_labels = String[]
            gene_expr_quantiles = Int[]
            quantile_legend = String[]
        end
    end

    if isempty(gene_expr_quantile_labels)

        println("ERROR: No expression quantiles could be created.")
        return
    end

    n_quantiles = length(unique(gene_expr_quantiles))

    expr_df_quantiles = @pipe copy(expr_df) |> insertcols!(_, :Quantile => gene_expr_quantiles) |> sort(_, :Quantile)

    for (r, sample_inds) in enumerate(sample_ind_groups)

        sample_name = gene_list[1].samples[sample_inds[1]]
        quantile_means = [Vector{Float64}() for i in 1:n_quantiles]

        for d in 1:n_quantiles

            gene_ids = expr_df_quantiles[expr_df_quantiles.Quantile .== d, :GeneID]
            gene_inds = findall([gene.id in gene_ids for gene in gene_list])

            for i in eachindex(gene_inds)

                # Get the average count within the specified region for each sample in sample_inds
                gene_sigs = [getsiginrange(gene_list[gene_inds[i]], region_list[r], sample_ind) for sample_ind in sample_inds]
                if !any(ismissing.(gene_sigs))
                    pos_count_gene = mean(reduce(vcat, gene_sigs))

                end

                if !any(ismissing.(gene_sigs))
                    push!(quantile_means[d], pos_count_gene / global_means[r])
                
                else
                    @warn "skipping $(gene_list[gene_inds[i]].id)
                    Does not have data in the specified region for samples $(gene_list[gene_inds[i]].samples[sample_inds])."
                
                end
            end
        end

        quant_ranges = String[]

        for d in 1:n_quantiles
            gene_inds = findall(gene_expr_quantiles .== d)
            quant_min = minimum(expr_df.Avg[gene_inds])
            quant_max = maximum(expr_df.Avg[gene_inds])
            
            if d == n_quantiles
            
                quant_range = @sprintf("[%.2f - %.2f]", quant_min, quant_max)
            else
            
                quant_range = @sprintf("[%.2f - %.2f)", quant_min, quant_max)
            end
            
            push!(quant_ranges, quant_range)
        end

        if box_plots
            vals = []
            quantile_nums = []
            for (i, mean_dist) in enumerate(quantile_means)
                if !isempty(mean_dist)
                    vals = vcat(vals, mean_dist)
                    quantile_nums = vcat(quantile_nums, fill(i, length(mean_dist)))
                end
            end
            # Create a box plot for the values across quantiles
            temp_df = DataFrame(
                :X => quantile_nums,
                :Y => vals
            )
            if horizontal
                bar_plt = plot(temp_df, x=:Y, y=:X, kind="box", orientation = "h")
            else
                bar_plt = plot(temp_df, x=:X, y=:Y, kind="box")
            end
        elseif isnothing(y_range)
            bar_plt = bar(
                x=[mean(mean_dist) for mean_dist in quantile_means],
                y=quant_ranges,
                orientation = horizontal ? "h" : "v"
            )
        elseif length(y_range) == 2
            if horizontal
                bar_plt = bar(
                    x=[mean(mean_dist) for mean_dist in quantile_means],
                    y=quant_ranges,
                    orientation = "h",
                    xaxis=attr(range=y_range)
                )
            else
                bar_plt = bar(
                    x=[mean(mean_dist) for mean_dist in quantile_means],
                    y=quant_ranges,
                    orientation = "v",
                    yaxis=attr(range=y_range)
                )
            end
        else
            @warn "Invalid y_range. Must be either nothing or a vector of length 2."
            bar_plt = bar(
                x=[mean(mean_dist) for mean_dist in quantile_means],
                y=quant_ranges,
                orientation = horizontal ? "h" : "v"
            )
        end

        if save_plots
            savefig(box_plots ? bar_plt : plot(bar_plt), joinpath(plot_save_dir, "$(sample_name)_average_enrichment_vs_expr_$(n_quantiles)-quantile_bar.html"))

        end

        if return_figs
            push!(fig_vec, bar_plt)

        else
            display(plot(bar_plt))
        
        end

        push!(p_values, KruskalWallisTest(quantile_means...))
        push!(means_vecs, quantile_means)
    end

    if return_figs
        return fig_vec, p_values, means_vecs

    end

    return p_values, means_vecs
end

function bootstrapavgenrich(gene_list::Vector{F}, sample_ind::Int, sample_size::Int, n_bootstraps::Int; upstream_lim::Int=1999, downstream_lim::Int=2000) where F <: Feature

    region_length = length(-upstream_lim:downstream_lim)
    bootstrap_mat = zeros(region_length, n_bootstraps)

    for i in 1:n_bootstraps

        bootstrap_mat[:,i] = getavgenrich(gene_list[rand(1:length(gene_list), sample_size)], sample_ind, upstream_lim=upstream_lim, downstream_lim=downstream_lim)
    end

    return bootstrap_mat
end

function getpval(test_group::Matrix{Int}, theta0_mat::Matrix{Float64}; adjust_p_value::Bool=true, p_val_mag::Bool=false) 

    p_vals = Float64[]
    
    if p_val_mag
        
        signs = Float64[]

        for i in 1:size(test_group)[1]

            test_result = ChisqTest(test_group[i,:], theta0_mat[i,:])
            push!(signs, sign(test_result.thetahat[1] - test_result.theta0[1]))
            temp_p_val = pvalue(test_result)
            push!(p_vals, temp_p_val)
        end

        if adjust_p_value

            p_vals = MultipleTesting.adjust(p_vals, BenjaminiHochberg())
        end

        return (1 .- p_vals) .* signs
    else

        for i in 1:size(test_group)[1]

            temp_p_val = pvalue(ChisqTest(test_group[i,:], theta0_mat[i,:]))
            push!(p_vals, temp_p_val)
        end

        if adjust_p_value

            p_vals = MultipleTesting.adjust(p_vals, BenjaminiHochberg())
        end

        return p_vals
    end
end

function p_val_enrich(feat_list::Vector{F}, null_list::Vector{G}, sample_inds::Union{Int, Vector{Int}};
                        four_regions::Bool=false, 
                        region_lims::Vector{R}=[1:1500,(1501:2100),(2101:2600),(2601:4000)]) where {F <: Feature, G <: Feature, R <: UnitRange}

    # Average the count across replicates/related samples:
    enrich_counts = Int.(round.(mean([getcountenrich(feat_list, i) for i in [sample_inds...]])))
    theta0 = mean([getavgenrich(null_list, i) for i in [sample_inds...]]) # theta0

    if four_regions
    
        enrich_counts_1 = Int.(round.(map(v -> mean(v), [enrich_counts[first(lim):last(lim),1] for lim in region_lims])))
        enrich_counts_2 = Int.(round.(map(v -> mean(v), [enrich_counts[first(lim):last(lim),2] for lim in region_lims])))
        enrich_counts = [enrich_counts_1 enrich_counts_2]
        theta0 = map(v -> mean(v), [theta0[first(lim):last(lim)] for lim in region_lims])
    end

    theta0 = [theta0 (1 .- theta0)]
    # return (enrich_counts, theta0) # TODO: REMOVE

    p_val_mat = getpval(enrich_counts, theta0)
    return p_val_mat
end

function averagecoveragebin(feature_list::Vector{F}, sample_ind::Int; region_cov::Bool=false) where {F<:Feature}

    feature_cov_avg = zeros(length(feature_list))

    if region_cov

        # TODO
        error("under construction")
    else

        for (i, feature) in enumerate(feature_list)

            temp_avg = sum(feature.binsignals[sample_ind]) / length(feature.binsignals[sample_ind])
            feature_cov_avg[i] = temp_avg
        end
    end

    return feature_cov_avg
end

function symmetricenrich(repeat_list::Vector{Repeat}, signal_ind::Int; peak_data::Bool=true, region_length::Int=1000)

    if isodd(region_length)

        upstream_limit = Int(floor(region_length/2))
        downstream_limit = upstream_limit
    else

        upstream_limit = Int(region_length/2)
        downstream_limit = upstream_limit - 1
    end

    n_counted = 0
    enrich_vec = zeros(region_length)

    for repeat_elem in repeat_list

        repeat_length = repeat_elem.repeat_end - repeat_elem.repeat_start
        
        if repeat_length >= region_length
        
            if isodd(repeat_length)

                center_coord = Int(ceil(repeat_length/2))
            else

                center_coord = Int(repeat_length/2)
            end

            downstream_limit_temp = center_coord + downstream_limit
            upstream_limit_temp = (center_coord - upstream_limit)
            signal_range = upstream_limit_temp:downstream_limit_temp

            forward_signal = repeat_elem.binsignals[signal_ind][signal_range]
            reverse_signal = reverse(forward_signal)
            average_signal = (forward_signal .+ reverse_signal) ./ 2
            enrich_vec .+= average_signal
            n_counted += 1
        end
    end

    return (enrich_vec ./ n_counted, n_counted)
end

function symmetricenrich(gene_list::Vector{Gene}, signal_ind::Int; peak_data::Bool=true, region_length::Int=1000)

    if isodd(region_length)

        upstream_limit = Int(floor(region_length/2))
        downstream_limit = upstream_limit
    else

        downstream_limit = Int(region_length/2)
        upstream_limit = downstream_limit - 1
    end

    n_counted = 0
    enrich_vec = zeros(region_length)

    for gene in gene_list

        gene_length = gene.gene_end - gene.gene_start
        
        if gene_length >= region_length
        
            if isodd(gene_length)

                center_coord = Int(ceil(gene_length/2))
            else

                center_coord = Int(gene_length/2)
            end

            downstream_limit_temp = center_coord + downstream_limit
            upstream_limit_temp = (center_coord - upstream_limit)
            signal_range = upstream_limit_temp:downstream_limit_temp

            forward_signal = gene.binsignals[signal_ind][signal_range]
            reverse_signal = reverse(forward_signal)
            average_signal = (forward_signal .+ reverse_signal) ./ 2
            enrich_vec .+= average_signal
            n_counted += 1
        end
    end

    return (enrich_vec ./ n_counted, n_counted)
end

function intergenic_dist(ref_genome)

    error("under construction")
end

"""Given vectors X and Y calculate ρ for N permutations of Y"""
function perm_cor_2side(X::Vector{Float64}, Y::Vector{Float64}, N::Int=10000)
    original_cor = cor(X, Y)
    Yperm = copy(Y)
    perm_cors = Float64[]
    for _ in 1:N
        shuffle!(Yperm)
        push!(perm_cors, abs(cor(X, Yperm)))
    end
    return (
        count(c -> c >= abs(original_cor), perm_cors) / (N + 1),
        original_cor
    )
end

"""function to get the correlation between enrichment in significant regions and dS"""
function get_cor(paralog_df::DataFrame, 
                 gene_range::GeneRange, 
                 sample_inds::Vector{Int}, 
                 genome::RefGenome, 
                 global_mean::Float64)
    XS = []
    YS = []
    for sample_ind in sample_inds
        genes = get(genome, collect(paralog_df.GeneID))
        paralogs = get(genome, collect(paralog_df.ParalogID))
        enrich_vals_gene = [!siginrange(gene, gene_range) ? missing : mean(getsiginrange(gene, gene_range, sample_ind)) for gene in genes]
        enrich_vals_paralog = [!siginrange(paralog, gene_range) ? missing : mean(getsiginrange(paralog, gene_range, sample_ind)) for paralog in paralogs]
        enrich_means = [(!ismissing(pair[1]) && !ismissing(pair[2])) ? mean(pair) / global_mean : missing for pair in zip(enrich_vals_gene, enrich_vals_paralog)]
        pairs_filtered = [pair for pair in zip(enrich_means, paralog_df.dS) if !ismissing(pair[1]) && !ismissing(pair[2])]
        xs = [pair[1] for pair in pairs_filtered]
        ys = [pair[2] for pair in pairs_filtered]
        push!(XS, xs)
        push!(YS, ys)
    end
    XS = reduce(hcat, XS)
    YS = reduce(hcat, YS)
    XS = vec(mean(XS, dims=2))
    YS = vec(mean(YS, dims=2))
    return perm_cor_2side(XS, YS)
end

"""function to get the correlation between enrichment in significant regions and dS"""
function get_cor_expr(expr_df::DataFrame, 
                      gene_range::GeneRange, 
                      sample_inds::Vector{Int}, 
                      genome::RefGenome)
    XS = []
    YS = []
    for sample_ind in sample_inds
        genes = get(genome, collect(expr_df.GeneID))
        xs = [getsiginrange(gene, gene_range, sample_ind) for gene in genes]
        ys = [expr_df.Avg[i] for (i, gene_sig) in enumerate(xs) if !ismissing(gene_sig)]
        push!(XS, mean.(collect(skipmissing(xs))))
        push!(YS, ys)
    end
    XS = reduce(hcat, XS)
    YS = reduce(hcat, YS)
    XS = vec(mean(XS, dims=2))
    YS = vec(mean(YS, dims=2))
    return perm_cor_2side(XS, YS)
end


# TESTS:
#=
test_gene1 = Gene(missing,
                  missing,
                  "test1",
                  "test1",
                  '+',
                  missing,
                  missing,
                  missing,
                  missing,
                  missing,
                  missing,
                  missing,
                  missing,
                  missing,
                  [],
                  missing,
                  1,
                  20,
                  nothing,
                  [],
                  [],
                  nothing
)  
=#