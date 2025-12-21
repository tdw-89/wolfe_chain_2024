if @isdefined ChromData
    
    nothing
else

    include("genomic_data.jl")
end

include("misc_utils.jl")

function addrepeats!(ref_genome::RefGenome, stk_file::String)

    temp_handle = open(stk_file, "r")

    temp_chrom = missing
    temp_start = missing
    temp_end = missing
    temp_type = missing
    temp_family = missing

    for line in eachline(temp_handle)
        
        if contains(line, "#=GF ID")

            temp_family = string(last(split(line, ' ')))
        elseif contains(line, "#=GF TP")
            
            temp_type = string(last(split(line, ' ')))
        elseif !(isnothing(match(r":[0-9]+-[0-9]+", line)))
            
            temp_chrom = string(split(line, ":")[1])
            range_str = match(r":[0-9]+-[0-9]+", line).match[2:end]
            range_split = string.(split(range_str, "-"))
            range_1 = parse(Int, range_split[1])
            range_2 = parse(Int, range_split[2])
            temp_start = min(range_1, range_2)
            temp_end = max(range_1, range_2)
            
            if haskey(ref_genome.scaffolds, temp_chrom)
            
                temp_chrom = ref_genome.scaffolds[temp_chrom]
            else

                ref_genome.scaffolds[temp_chrom] = Scaffold(temp_chrom, missing, Gene[], Repeat[], missing, missing, missing)
                temp_chrom = ref_genome.scaffolds[temp_chrom]
            end
            
            temp_repeat = Repeat(temp_chrom, 
                                missing,
                                Vector{Region}(),
                                temp_family, 
                                temp_type, 
                                temp_start, 
                                temp_end, 
                                missing,
                                Vector{Union{Vector{UInt8},Vector{UInt16},Vector{UInt32},Vector{UInt64}}}(),
                                Vector{BitVector}(),
                                Vector{String}())
            push!(ref_genome.repeats, temp_repeat)
            push!(temp_repeat.scaffold.repeats, temp_repeat)
        elseif line == "//"
            
            temp_chrom = missing
            temp_start = missing
            temp_end = missing
            temp_type = missing
            temp_family = missing
        end
    end

    close(temp_handle)
end

function move_to_repeats!(ref_genome::RefGenome, gene_list::Vector{Gene})

    for gene in gene_list

        temp_repeat = Repeat(
            gene.scaffold,
            missing,
            missing,
            missing,
            missing,
            gene.gene_start,
            gene.gene_end,
            missing,
            nothing,
            Vector{BitVector}(),
            Vector{String}()
        )

        push!(ref_genome.repeats, temp_repeat)
        push!(gene.scaffold.repeats, temp_repeat)
    end

    all_ids = [gene.id for gene in gene_list]
    filter!(gene -> gene.id ∉ all_ids, ref_genome.genes[2])
    filter!(gene -> gene ∉ all_ids, ref_genome.genes[1])
    for (scaffold_name, scaffold) in ref_genome.scaffolds

        filter!(gene -> gene.id ∉ all_ids, scaffold.genes)
    end
end

function findoverlappinggenes(repeat_elem::Repeat)

    repeat_start = repeat_elem.repeat_start
    repeat_end = repeat_elem.repeat_end
    overlap_ids = String[]
    repeat_type = Vector{Union{String,Missing}}()
    repeat_family = Vector{Union{String,Missing}}()

    for gene in repeat_elem.scaffold.genes

        gene_start = gene.gene_start
        gene_end = gene.gene_end

        if hasoverlap(repeat_start, gene_start, repeat_end, gene_end)

            push!(overlap_ids, gene.id)
            push!(repeat_type, repeat_elem.type)
            push!(repeat_family, repeat_elem.family)
        end
    end

    if !isempty(overlap_ids)

        return overlap_ids, repeat_type, repeat_family
    else

        return missing
    end
end

function parsestk(stk_file::String)

    temp_handle = open(stk_file, "r")

    temp_chrom = missing
    temp_start = missing
    temp_end = missing
    temp_type = missing
    temp_family = missing
    repeat_df = DataFrame(Chromosome=Vector{Union{String, Missing}}(), 
                            Start=Vector{Union{Int, Missing}}(), 
                            End=Vector{Union{Int, Missing}}(), 
                            Family=Vector{Union{String, Missing}}(), 
                            Type=Vector{Union{String, Missing}}())

    for line in eachline(temp_handle)
        
        if contains(line, "#=GF ID")

            temp_family = string(last(split(line, ' ')))
        elseif contains(line, "#=GF TP")
            
            temp_type = string(last(split(line, ' ')))
        elseif !(isnothing(match(r":[0-9]+-[0-9]+", line)))
            
            temp_chrom = string(split(line, ":")[1])
            range_str = match(r":[0-9]+-[0-9]+", line).match[2:end]
            range_split = string.(split(range_str, "-"))
            range_1 = parse(Int, range_split[1])
            range_2 = parse(Int, range_split[2])
            temp_start = min(range_1, range_2)
            temp_end = max(range_1, range_2)
            
            push!(repeat_df, (temp_chrom, temp_start, temp_end, temp_family, temp_type))
        elseif line == "//"
            
            temp_chrom = missing
            temp_start = missing
            temp_end = missing
            temp_type = missing
            temp_family = missing
        end
    end

    close(temp_handle)
    return repeat_df
end

function convert_to_repeats!(ref_genome::RefGenome, 
                             repeat_df::DataFrame;
                             allow_missing_scaffolds::Bool=false)

    for i in 1:nrow(repeat_df)

        if allow_missing_scaffolds

            if !haskey(ref_genome.scaffolds, repeat_df[i, :Chromosome])
                ref_genome.scaffolds[repeat_df[i, :Chromosome]] = Scaffold(repeat_df[i, :Chromosome], missing, Gene[], Repeat[], missing, missing, missing)
           
            end
        end
        
        temp_repeat = Repeat(
            ref_genome.scaffolds[repeat_df[i, :Chromosome]],
            missing,
            missing,
            repeat_df[i, :Family],
            repeat_df[i, :Type],
            repeat_df[i, :Start],
            repeat_df[i, :End],
            missing,
            nothing,
            Vector{BitVector}(),
            Vector{String}()
        )
        push!(temp_repeat.scaffold.repeats, temp_repeat)
        push!(ref_genome.repeats, temp_repeat)
        
    end
end

"""This function parses a Dfam '.embl' file and returns a DataFrame of the parsed data.
Info on the EMBL format: https://bibiserv.cebitec.uni-bielefeld.kw/sadr/data_formats/embl_df.html"""
function parseembl(embl_file::String)

    temp_handle = endswith(embl_file, ".gz") ? GzipDecompressorStream(open(embl_file, "r")) : open(embl_file, "r")
    
    ids, acs, kws = String[], String[], String[]

    id, ac, kw = missing, missing, missing
    for line in eachline(temp_handle)
        if startswith(line, "ID")
            id = split(line)[2][1:end-1] # removes trailing ';'
            ac = missing
            kw = missing
        elseif startswith(line, "AC")
            ac = split(line)[2][1:end-1] # removes trailing ';'
        elseif startswith(line, "KW")
            kw = join(split(line)[2:end])
        elseif line == "//"
            push!(ids, id)
            push!(acs, ac)
            push!(kws, kw)
            id, ac, kw = missing, missing, missing
        end
    end

    close(temp_handle)
    return DataFrame(ID=ids, AC=acs, KeyWords=kws)
end

"""This function parses a Dfam '.hits' file and returns a list of Repeat objects."""
function parsehits(hits_file::String, families_file::String; eval_cutoff::Float64=0.001) # <- NOTE: 0.001 is the cutoff used by default in dfam.org's sequence search tool

    println("Loading hits file")
    if endswith(hits_file, ".gz")
        hits_df = CSV.read(GzipDecompressorStream(open(hits_file, "r")), DataFrame)

    else
        hits_df = CSV.read(hits_file, DataFrame)

    end

    filter!(row -> row["e-value"] <= eval_cutoff, hits_df)

    # hits_df = hits_df[1:100000,:] # DEBUG REMOVE

    # Match the family accession numbers to the keywords from the families file
    families_df = parseembl(families_file)
    insertcols!(hits_df, 3, :KeyWords => "")
    println("Finding matching family keywords...                                              ↓")
    print("[")
    n_rows = nrow(hits_df)
    inc = n_rows ÷ 80
    inc = inc == 0 ? 1 : inc
    
    for i in 1:n_rows
        temp_id = hits_df[i, :family_acc]
        temp_kw = families_df[findfirst(id -> contains(temp_id, id), families_df.ID), :KeyWords]
        hits_df[i, :KeyWords] = temp_kw
        if i % inc == 0
            print("=")
        end
    end
    print("]\n")

    return hits_df
end