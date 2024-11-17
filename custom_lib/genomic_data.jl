# DEPENDENCIES #
using StatsBase
using DataFrames
using XAM
using CSV
using CodecZlib

if !@isdefined RefGenome
    include("genome_types.jl")

end

if !@isdefined PeakData
    include("./narrowPeak_parsing/PeakData.jl")

end

# TYPES #

mutable struct ChromData
    name::String
    signal::Union{Vector{UInt16}, Vector{Float64}, BitVector}

end

mutable struct SampleData
    name::String
    chroms::Vector{ChromData}
    n_reads::R where R <: Real

end

mutable struct Experiment    
    name::String
    samples::Vector{SampleData}

end

mutable struct PeakWarnings
    counts::Vector{Int}
    switches::Vector{Bool}

end

PeakWarnings() = PeakWarnings([0, 0], [true, true])

# FUNCTIONS #

""" This function takes a (full) path to a bam file, reads each record in the bam file,\n
and modifies the provided dictionary with chromosome names mapped to integer vectors,\n 
which contain the counts for # of reads overlapping each nucleotide. The vector indices\n 
correspond to the position on the respective chromosome."""
function getallrecords(bam_file_path::String; mapq_signal::Bool=false)

    bam_reader = open(BAM.Reader, bam_file_path)
    chroms = bam_reader.refseqnames
    chrom_lengths = bam_reader.refseqlens
    count_vec_list = Vector{ChromData}(undef, length(chroms))
    total_read_count = 0

    for i in eachindex(chroms)

        if mapq_signal
        
            count_vec_list[i] = ChromData(chroms[i], zeros(chrom_lengths[i]))
        else
        
            count_vec_list[i] = ChromData(chroms[i], UInt16.(zeros(chrom_lengths[i])))
        end
    end
    
    temp_record = BAM.Record()

    while !eof(bam_reader)

        empty!(temp_record)
        read!(bam_reader, temp_record)
        total_read_count += 1
        
        if BAM.ismapped(temp_record)

            start_pos = BAM.position(temp_record)
            end_pos = BAM.rightposition(temp_record)
            chrom_name = BAM.refname(temp_record)

            ind = findfirst([chrom_name == count_vec_list[j].name for j in 1:length(count_vec_list)])

            if mapq_signal

                mapq = mean(BAM.mappingquality(temp_record))
                
                @inbounds for k in start_pos:end_pos

                    count_vec_list[ind].signal[k] += mapq
                end
            else

                @inbounds for k in start_pos:end_pos

                    count_vec_list[ind].signal[k] += 1
                end
            end

        else
            nothing # TO DO: print/store name for unmapped records?
        end
    end
    
    close(bam_reader)

    return count_vec_list, total_read_count
end

"""
This function takes a vector of bam file paths, calls 'getallrecords' for\n
each one, and returns an 'Experiment' object containing the vectors of read\n
counts for each chromosome in each bam file.\n
"""
function getallcountvectors(bam_file_list::Vector{String}; exp_name::String="seq_exp", mapq_signal::Bool=false)

    sample_names = basename.(bam_file_list)
    bam_count_list = Experiment(exp_name, Vector{SampleData}(undef, length(bam_file_list)))

    for i in 1:length(bam_file_list)

        ret_vec, total_counts = getallrecords(bam_file_list[i], mapq_signal=mapq_signal)
        bam_count_list.samples[i] = SampleData(sample_names[i], 
                                               Vector{ChromData}(undef, length(ret_vec)),
                                               total_counts)

        for k in eachindex(ret_vec)

            bam_count_list.samples[i].chroms[k] = ret_vec[k]
        end
    end

    return bam_count_list
end

""" 
    Given an 'Experiment' object and a list of replicates, this function will return another 
    'Experiment' object with the average signal for each replicate group.\n
    The list can be an n column file (i.e. "replicate 1"\t"replicate 2"\t"replicate 3"...\t"replicate n")
    with no header or a DataFrame with n columns.
"""
function average_bam_replicate_groups(bam_data::Experiment; replicate_list::Union{String, DataFrame, Nothing}=nothing)
    n_chroms = length(bam_data.samples[1].chroms)
    @assert all([length(bam_data.samples[i].chroms) == n_chroms for i in 2:length(bam_data.samples)])

    if isa(replicate_list, String)
        replicate_list = CSV.read(replicate_list, DataFrame; header=false)

    elseif isnothing(replicate_list)
        sample_names = [bam_data.samples[i].name for i in eachindex(bam_data.samples)]
        replicate_list = DataFrame(["$i" => sample_name for (i, sample_name) in enumerate(sample_names)])

    end

    avg_samples = SampleData[]

    for i in 1:nrow(replicate_list)
        replicate_names = replicate_list[i,:]
        replicate_vec = [rep_sample for rep_sample in bam_data.samples if rep_sample.name in replicate_names]
        push!(avg_samples, average_bam_replicates(replicate_vec))

    end

    return Experiment("avg_" * bam_data.name, avg_samples)
end

function average_peak_replicate_groups(peak_data::Experiment; replicate_list::Union{String, DataFrame, Nothing}=nothing)
    if isa(replicate_list, String)
        replicate_list = CSV.read(replicate_list, DataFrame; header=false)

    elseif isnothing(replicate_list)
        sample_names = [peak_data.samples[i].name for i in eachindex(peak_data.samples)]
        replicate_list = DataFrame(["$i" => sample_name for (i, sample_name) in enumerate(sample_names)])

    end

    avg_samples = SampleData[]

    for i in 1:nrow(replicate_list)
        replicate_names = replicate_list[i,:]
        replicate_vec = [rep_sample for rep_sample in peak_data.samples if rep_sample.name in replicate_names]
        push!(avg_samples, average_peak_replicates(replicate_vec))

    end

    return Experiment("avg_" * peak_data.name, avg_samples)

end

function average_bam_replicates(replicate_vec::Vector{SampleData})
    # Get total count ratios:
    count_mean = mean([replicate_vec[i].n_reads for i in eachindex(replicate_vec)])
    count_ratios = [count_mean / replicate_vec[i].n_reads for i in eachindex(replicate_vec)]
    return_chroms = [ChromData(replicate_vec[1].chroms[i].name, zeros(length(replicate_vec[1].chroms[i].signal))) for i in eachindex(replicate_vec[1].chroms)]

    for chrom_n in eachindex(return_chroms)
        for rep_n in eachindex(replicate_vec)
            return_chroms[chrom_n].signal .+= replicate_vec[rep_n].chroms[chrom_n].signal .* count_ratios[rep_n]

        end

        return_chroms[chrom_n].signal ./= length(replicate_vec)

    end

    return SampleData("avg_" * replicate_vec[1].name, return_chroms, count_mean)

end

function average_peak_replicates(replicate_vec::Vector{SampleData})
    n_chroms = length(replicate_vec[1].chroms)
    @assert all([length(replicate_vec[i].chroms) == n_chroms for i in 2:length(replicate_vec)]) "All samples must have the same number of chromosomes"
    @assert all([all([length(replicate_vec[1].chroms[j].signal) == length(replicate_vec[i].chroms[j].signal) for j in 1:n_chroms]) for i in 2:length(replicate_vec)]) "All samples must have the same chromosome lengths"
    return_chroms = [ChromData(replicate_vec[1].chroms[i].name, zeros(length(replicate_vec[1].chroms[i].signal))) for i in eachindex(replicate_vec[1].chroms)]

    for chrom_n in eachindex(return_chroms)

        for rep_n in eachindex(replicate_vec)
            return_chroms[chrom_n].signal .+= replicate_vec[rep_n].chroms[chrom_n].signal

        end

        return_chroms[chrom_n].signal ./= length(replicate_vec)

    end

    return SampleData("avg_" * replicate_vec[1].name, return_chroms, 1)

end

function addpeak!(chrom_dict::Dict{String, BitVector}, chrom::AbstractString, peak_start::Int, peak_end::Int, warnings::PeakWarnings)

    if haskey(chrom_dict, chrom)
        try
            chrom_dict[chrom][peak_start:peak_end] .= 1

        catch BoundsError
            if warnings.switches[1]
                @warn "Bounds error for peak on chrom $(chrom)"
                warnings.counts[1] += 1

                if warnings.counts[1] ==10
                    @warn "Too many bounds error warnings were issued. Suppressing further warnings..."
                    warnings.switches[1] = false

                end
            end

            chrom_dict[chrom][clamp(peak_start, 1, length(chrom_dict[chrom])):clamp(peak_end, 1, length(chrom_dict[chrom]))] .= 1

        end

    elseif haskey(chrom_dict, chrom[4:end])

        if warnings.switches[2]
            @warn "Chromosome '$(chrom)' not found in chromosome lengths file, using '$(chrom[4:end])' instead."
            warnings.counts[2] += 1

            if warnings.counts[2] ==10
                @warn "Too many chromosome name warnings were issued. Suppressing further warnings..."
                warnings.switches[2] = false

            end
        end
        
        try 
            chrom_dict[chrom[4:end]][peak_start:peak_end] .= 1

        catch BoundsError
            if warnings.switches[1]
                @warn "Bounds error for peak on chrom $(chrom)"
                warnings.counts[1] += 1

                if warnings.counts[1] ==10
                    @warn "Too many bounds error warnings were issued. Suppressing further warnings..."
                    warnings.switches[1] = false

                end
            end

            chrom_dict[chrom[4:end]][clamp(peak_start, 1, length(chrom_dict[chrom[4:end]])):clamp(peak_end, 1, length(chrom_dict[chrom[4:end]]))] .= 1

        end
    else
        error("Chromosome $(chrom) not found in chromosome lengths file")

    end
end

"""Load narrowPeak file data as binary peak vectors."""
function binpeaks(peak_files::Union{String, Vector{String}}, chrom_lengths_file::Union{String, Nothing}=nothing)

    if typeof(peak_files) == String

        peak_files = [peak_files]
    end

    chrom_lengths_df = CSV.read(chrom_lengths_file, DataFrame, header=false)
    chrom_lengths_df[!, 1] = string.(chrom_lengths_df[!, 1])
    rename!(chrom_lengths_df, ["chrom","length"])

    ret_experiment = Experiment("peak_data", Vector{SampleData}(undef, length(peak_files)))
    
    chrom_dict = Dict{String, BitVector}()
    warnings = PeakWarnings()
    
    for i in eachindex(peak_files)
        
        for j in 1:nrow(chrom_lengths_df)

            chrom_dict[chrom_lengths_df.chrom[j]] = BitVector(zeros(chrom_lengths_df.length[j]))
        end

        temp_df = DataFrame()

        try
            temp_df = PeakData.LoadNarrowPeak(peak_files[i])
        catch ArgumentError
        
            colnames = [:chrom,
                        :start,
                        :end,
                        :peakName,
                        :score,
                        :strand,
                        :signalValue,
                        :pValue,
                        :qValue,
                        :peakPoint]
            
            temp_df = endswith(peak_files[i], ".gz") ? 
                        CSV.read(GzipDecompressorStream(open(peak_files[i], "r")), DataFrame, header=false) : 
                        CSV.read(peak_files[i], DataFrame, header=false)
            
            rename!(temp_df, colnames)
        end

        for peak in eachrow(temp_df)

            addpeak!(chrom_dict, peak.chrom, peak.start + 1, peak.end, warnings) # NOTE: .narrowPeak (and .bed files in general) use 0-based half-open intervals
        end

        temp_name = String(split(basename(peak_files[i]), ".")[1])
        # temp_name = reduce(*, temp_name[1:length(temp_name) - 1])
        # temp_name = split(temp_name, "_")[1:2]
        # temp_name = temp_name[1] * "_" * temp_name[2]
        temp_sample = SampleData(temp_name, ChromData[], -1)

        for (chrom_name, signal_data) in chrom_dict

            push!(temp_sample.chroms, ChromData(chrom_name, signal_data))
        end

        ret_experiment.samples[i] = temp_sample
    end

    return ret_experiment
end

function binpeaks(narrow_peak_dir::String, chrom_lengths_file::Union{String, Nothing}=nothing)

    peak_files = readdir(narrow_peak_dir, join=true)
    peak_files = peak_files[contains.(peak_files, r".narrowPeak$") .|| contains.(peak_files, r".bed$")]
    return binpeaks(peak_files, chrom_lengths_file)
end

# NOTE: This function should be calling a HOMER specific peak file parser, to be located in the PeakData module.
# Currently, it just cheats and uses CSV.jl and some hard coded column names.
function binpeakshomer(peak_files::Vector{String}, chrom_lengths_file::Union{String, Nothing}=nothing; gff_files::Union{String, Nothing, Vector{String}}=nothing)

    chrom_lengths_df = DataFrame()

    if isnothing(chrom_lengths_file)

        if isnothing(gff_files)

            error("Must provide either a file containing chromosome lengths or a GFF file.") # COPILOT: Error message was written by github copilot after entering 'er'
        else

            for gff_file in gff_files

                chrom_lengths_df = vcat(chrom_lengths_df, getchromlengths(gff_file))
            end

            rename!(chrom_lengths_df, ["chrom","length"])
        end
    else

        chrom_lengths_df = CSV.read(chrom_lengths_file, DataFrame, header=false)
        chrom_lengths_df[!, 1] = string.(chrom_lengths_df[!, 1])
        rename!(chrom_lengths_df, ["chrom","length"])
    end

    ret_experiment = Experiment("peak_data", Vector{SampleData}(undef, length(peak_files)))
    chrom_dict = Dict{String, BitVector}()
    warnings = PeakWarnings()

    for i in 1:length(peak_files)
        
        for j in 1:nrow(chrom_lengths_df)

            chrom_dict[chrom_lengths_df.chrom[j]] = BitVector(zeros(chrom_lengths_df.length[j]))
        end

        temp_df = CSV.read(peak_files[i], header=7, DataFrame)
        filter!(row -> row."p-value vs Control" <= 0.01, temp_df)
        temp_df = temp_df[:, ["#PeakID", "chr", "start", "end"]]
        rename!(temp_df, ["peakName", "chrom", "start", "end"])

        for peak in eachrow(temp_df)
        
            addpeak!(chrom_dict, String(peak.chrom), peak.start + 1, peak.end, warnings) # NOTE: .narrowPeak (and .bed files in general) use 0-based half-open intervals
            
        end

        temp_name = split(basename(peak_files[i]), ".")[1]
        temp_sample = SampleData(temp_name, ChromData[], -1)

        for (chrom_name, signal_data) in chrom_dict

            push!(temp_sample.chroms, ChromData(chrom_name, signal_data))
        end

        ret_experiment.samples[i] = temp_sample
    end

    return ret_experiment
end


"""Add peak data to a genome object."""
function addtogenes!(genome::RefGenome, experiment::Experiment; peak_data::Bool=true, regions::Bool=true)

    if peak_data
    
        if typeof(experiment.samples[1].chroms[1].signal) != BitVector

            error("If 'peak_data' is set to true, then the signal must be a BitVector")
        end
    end
    
    for gene in genome.genes[2]

        gene_start = gene.gene_start
        gene_end = gene.gene_end
        
        if ismissing(gene.scaffold)
            id = gene.id
            error("'Gene' $id is missing chromosome information (a 'Scaffold' object)")
        end

        chrom_name = gene.scaffold.name

        for sample_data in experiment.samples
            
            for chrom_data in sample_data.chroms
            
                if chrom_data.name == chrom_name
                
                    if peak_data

                        
                        signal_vec = chrom_data.signal[gene_start:gene_end]
                        push!(gene.binsignals, signal_vec)

                        if regions

                            if isempty(gene.regions)

                                # TO DO: add default region
                                error("adding regions inside 'addtogenes!' not implemented yet") # DEBUG REMOVE
                            else

                                region_start = gene.regions[1].region_start
                                region_end = gene.regions[1].region_end
                                region_signal_vec = chrom_data.signal[region_start:region_end]
                                push!(gene.regions[1].binsignals, region_signal_vec)
                            end
                        end

                        break
                    else
                        if typeof(chrom_data.signal[1]) == Float64

                            signal_vec = chrom_data.signal[gene_start:gene_end]
                        else
                            signal_vec = try UInt8.(chrom_data.signal[gene_start:gene_end])

                            catch InexactError

                                signal_vec = try UInt16.(chrom_data.signal[gene_start:gene_end])

                                catch InexactError

                                    signal_vec = try UInt32.(chrom_data.signal[gene_start:gene_end])

                                    catch InexactError

                                        signal_vec = UInt64.(chrom_data.signal[gene_start:gene_end])
                                    end
                                end
                            end
                        end

                        push!(gene.signals, signal_vec)
                        
                        if regions
                        
                            if isempty(gene.regions)

                                # TO DO: add default region
                                error("adding regions inside 'addtogenes!' not implemented yet") # DEBUG REMOVE
                            else

                                region_start = gene.regions[1].region_start
                                region_end = gene.regions[1].region_end

                                if typeof(chrom_data.signal[1]) == Float64

                                    region_signal_vec = chrom_data.signal[region_start:region_end]
                                else 

                                    region_signal_vec = try UInt8.(chrom_data.signal[region_start:region_end])

                                    catch InexactError
            
                                        region_signal_vec = try UInt16.(chrom_data.signal[region_start:region_end])
            
                                        catch InexactError
            
                                            region_signal_vec = try UInt32.(chrom_data.signal[region_start:region_end])
            
                                            catch InexactError
            
                                                region_signal_vec = UInt64.(chrom_data.signal[region_start:region_end])
                                            end
                                        end
                                    end
                                end

                                push!(gene.regions[1].signals, region_signal_vec)
                            end
                        end
                        
                        break
                    end
                end
            end

            sample_name = sample_data.name
            push!(gene.samples, sample_name)

            if regions

                push!(gene.regions[1].samples, sample_name)
            end
        end
    end
end

"""Add expression data to a genome object.\n
NOTE: this function assumes that the first column is the gene id column, all others are exprsesion values"""
function addexpression!(ref_genome::RefGenome, expr_data::DataFrame; total_expr::Bool=true)

    # Add expression data as 'RNA' objects to the genome object. Expression values
    # from different samples are added in the order of their respective columns.
    for gene in eachrow(expr_data)
        id = gene.GeneID
        gene_ind = findfirst(id .== ref_genome.genes[1])

        if !(isnothing(gene_ind))
            new_rna = RNA("all_transcripts", "mRNA", ref_genome.genes[2][gene_ind], Exon[], Intron[], names(gene[3:end]), missing, missing, Float64.(Array(gene[2:end])))
            pushfirst!(ref_genome.genes[2][gene_ind].rnas, new_rna)

        end
    end
end

function addtorepeats!(genome::RefGenome, experiment::Experiment; peak_data::Bool=false, add_to_region::Bool=false)

    if peak_data
    
        if typeof(experiment.samples[1].chroms[1].signal) != BitVector

            error("If 'peak_data' is set to true, then the signal must be a BitVector")
        end
    else

        if typeof(experiment.samples[1].chroms[1].signal) != Vector{UInt16}

            error("If 'peak_data' is set to false (default), then the signal must be a UInt16 vector")
        end
    end
    
    for repeat_elem in genome.repeats
        
        if ismissing(repeat_elem.scaffold) # TO DO: add 'repeat' to error message
            
            error("")
        end

        chrom_name = repeat_elem.scaffold.name
        repeat_start = repeat_elem.repeat_start
        repeat_end = repeat_elem.repeat_end

        for sample_data in experiment.samples
            
            for chrom_data in sample_data.chroms
            
                if chrom_data.name == chrom_name
                
                    if peak_data

                        if add_to_region

                            if isempty(repeat_elem.regions)

                                error("No method for adding default region to repeat element") # DEBUG REMOVE
                            else

                                region_start = repeat_elem.regions[1].region_start
                                region_end = repeat_elem.regions[1].region_end
                                region_signal_vec = chrom_data.signal[region_start:region_end]
                                push!(repeat_elem.regions[1].binsignals, region_signal_vec)
                            end
                        end
                        
                        repeat_signal_vec = chrom_data.signal[repeat_start:repeat_end]
                        push!(repeat_elem.binsignals, repeat_signal_vec)
                        break
                    else
                        

                        if isempty(repeat_elem.regions)

                            error("No method for adding default region to repeat element") # DEBUG REMOVE
                        end

                        if add_to_region

                            region_start = repeat_elem.regions[1].region_start
                            region_end = repeat_elem.regions[1].region_end

                            region_signal_vec = try UInt8.(chrom_data.signal[region_start:region_end])

                            catch InexactError

                                region_signal_vec = try UInt16.(chrom_data.signal[region_start:region_end])

                                catch InexactError

                                    region_signal_vec = try UInt32.(chrom_data.signal[region_start:region_end])

                                    catch InexactError

                                        region_signal_vec = UInt64.(chrom_data.signal[region_start:region_end])
                                    end
                                end
                            end

                            push!(repeat_elem.regions[1].signals, region_signal_vec)
                        end
                        
                        repeat_signal_vec = try UInt8.(chrom_data.signal[repeat_start:repeat_end])

                        catch InexactError

                            repeat_signal_vec = try UInt16.(chrom_data.signal[repeat_start:repeat_end])

                            catch InexactError

                                repeat_signal_vec = try UInt32.(chrom_data.signal[repeat_start:repeat_end])

                                catch InexactError

                                    repeat_signal_vec = UInt64.(chrom_data.signal[repeat_start:repeat_end])
                                end
                            end
                        end

                        push!(repeat_elem.signals, repeat_signal_vec)
                        break
                    end
                end
            end

            sample_name = sample_data.name
            push!(repeat_elem.samples, sample_name)

            if add_to_region
                
                push!(repeat_elem.regions[1].samples, sample_name)
            end
        end
    end
end

function has_expr(gene::Gene)

    if isempty(gene.rnas)

        return false
    end

    if isempty(gene.rnas[1].expression)

        return false
    end

    return true
end

function sample_to_bed_df(sample::SampleData)
    chroms = String[]
    starts = Int[]
    ends = Int[]
    scores = Int[]

    for chrom in sample.chroms
        rg_val_pairs = contiguous_values(chrom.signal)
        push!(chroms, [chrom.name for _ in rg_val_pairs]...)
        push!(starts, [first(first(rg_val_pair)) for rg_val_pair in rg_val_pairs]...)
        push!(ends, [last(first(rg_val_pair)) for rg_val_pair in rg_val_pairs]...)
        push!(scores, [round(Int, last(rg_val_pair) * 1000) for rg_val_pair in rg_val_pairs]...)
    
    end

    df = DataFrame("Chrom" => chroms, "Start" => starts, "End" => ends, "Score" => scores) 
    insertcols!(df, 4, :Name => string.(1:nrow(df)))
    return df

end

"""Given a Float64 vector, return a list of all contiguous stretches 
of zeros (as integer ranges)""" 
function contiguous_values(signal::Vector{Float64}, val::Real, val_complement::Bool=false)

    contig_ranges = Vector{UnitRange{Int}}()
    start_ind = 0
    end_ind = 0
    in_contig = false

    for i in eachindex(signal)

        if val_complement
            if signal[i] != val

                if !in_contig

                    start_ind = i
                    in_contig = true
                end
            else

                if in_contig

                    end_ind = i - 1
                    push!(contig_ranges, start_ind:end_ind)
                    in_contig = false
                end
            end
        else
            if signal[i] == val

                if !in_contig

                    start_ind = i
                    in_contig = true
                end
            else

                if in_contig

                    end_ind = i - 1
                    push!(contig_ranges, start_ind:end_ind)
                    in_contig = false
                end
            end
        end
    end

    return contig_ranges
end

"""Given a Float64 vector, return a list of all contiguous stretches of any non-zero value,
along with the value itself"""
function contiguous_values(signal::Vector{Float64})

    contig_ranges = Vector{Tuple{UnitRange{Int}, Float64}}()
    start_ind = 0
    end_ind = 0
    in_contig = false
    val = -1.0

    for i in eachindex(signal)
        if signal[i] != 0.0
            if !in_contig
                start_ind = i
                val = signal[i]
                in_contig = true
                
            elseif signal[i] != val
                push!(contig_ranges, (start_ind:i-1, val))
                start_ind = i
                val = signal[i]
            end
        else
            if in_contig
                push!(contig_ranges, (start_ind:i-1, val))
                in_contig = false
            end
        end
    end

    if in_contig
        push!(contig_ranges, (start_ind:length(signal), val))

    end

    return contig_ranges
end

# SHOW OVERLOADING #

Base.show(io::IO, x::ChromData) = begin

    chrom_name = x.name
    sig_len = length(x.signal)
    print(io, "Chromosome $chrom_name with a $sig_len bp signal")
end

Base.show(io::IO, x::Experiment) = begin

    exp_name = x.name
    n_samples = length(x.samples)
    print(io, "Experiment '$exp_name' with $n_samples samples")
end

Base.show(io::IO, x::SampleData) = begin

    sample_name = x.name
    n_chroms = length(x.chroms)
    print(io, "Sample '$sample_name' with $n_chroms sequences" * "$(x.n_reads == -1 ? "" : " and $(x.n_reads) reads")")
end