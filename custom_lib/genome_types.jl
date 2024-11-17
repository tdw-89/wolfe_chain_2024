using BioSequences

abstract type 

    Feature
end

abstract type

    RegElement <: Feature
end

# mutable struct Genome <: Feature

#     version::Union{VersionNumber, Int, Float64, Missing, String}
#     organism::Union{String,Missing}
#     scaffolds::Union{Vector{Feature},Missing}
# end

mutable struct Scaffold <: Feature

    name::Union{String, Char}
    contigs::Union{Vector{Feature}, Missing}
    genes::Union{Vector{Feature}, Missing}
    repeats::Union{Vector{Feature}, Missing}
    scaffold_start::Union{Int, Missing}
    scaffold_end::Union{Int, Missing}
    level::Union{String, Missing}
end

mutable struct Contig <: Feature

    scaffold::Union{Scaffold, Missing} # Parent Node
    junk::Union{Vector{Feature}, Missing}
    enhancers::Union{Vector{Feature}, Missing}
    genes::Union{Vector{Feature}, Missing}
    contig_start::Union{Int, Missing}
    contig_end::Union{Int, Missing}
end

struct Enhancer <: RegElement

    contig::Union{Contig, Missing} # Parent Node
end

struct Gene <: Feature 

    # Intrinsic features:
    scaffold::Union{Scaffold, String, Missing}
    contig::Union{Contig, String, Missing} # Parent Node
    name::Union{String, Missing}
    id::Union{String, Missing}
    strand::Union{Char, Missing}
    cres::Union{Vector{RegElement}, Missing}
    tss::Union{Int, Missing}
    tes::Union{Int, Missing}
    introns::Union{Vector{Feature}, Missing}
    exons::Union{Vector{Feature}, Missing}
    cds_start::Union{Int, Missing}
    cds_end::Union{Int, Missing}
    rnas::Union{Vector{Feature}, Missing}
    segments::Union{Vector{Feature}, Missing}
    regions::Union{Vector{Feature}, Missing}
    annotations::Union{Dict{String, Vector{Feature}}, Missing}
    gene_start::Union{Int, Missing}
    gene_end::Union{Int, Missing} 
    sequence::Union{B, Nothing} where B <: BioSequences.BioSequence

    # Experimental data:
    signals::Union{Vector{Union{Vector{UInt8},Vector{UInt16},Vector{UInt32},Vector{UInt64}, Vector{Float64}}}, Nothing}
    binsignals::Union{Vector{BitVector}, Nothing}
    samples::Union{Vector{String}, Nothing}
end 

struct Promoter <: RegElement 

    # Intrinsic features:
    genes::Union{Vector{Gene}, Missing}
    strand::Union{Char, Missing}
    segments::Union{Vector{Feature}, Missing}
    promoter_start::Union{Int, Missing}
    promoter_end::Union{Int, Missing}
    sequence::Union{B, Nothing} where B <: BioSequences.BioSequence

    # Experimental data:
    signals::Union{Vector{Union{Vector{UInt8},Vector{UInt16},Vector{UInt32},Vector{UInt64}, Vector{Float64}}}, Nothing}
    binsignals::Union{Vector{BitVector}, Nothing}
    samples::Union{Vector{String}, Nothing}
end

struct Intron <: Feature

    gene::Union{Gene, Missing} # Parent Node
    rna::Union{Feature, Missing} # Parent Node
    intron_start::Union{Int, Missing}
    intron_end::Union{Int, Missing}
end

struct Exon <: Feature

    gene::Union{Gene, Missing} # Parent Node
    rna::Union{Feature, Missing} # Parent Node
    exon_start::Union{Int, Missing}
    exon_end::Union{Int, Missing}
end

struct RNA <: Feature

    # Intrinsic features:
    id::Union{String, Missing}
    type::Union{String, Missing}
    gene::Union{Gene, Missing}
    exons::Union{Vector{Exon}, Missing}
    introns::Union{Vector{Intron}, Missing}
    samples::Union{Vector{String}, Missing}
    rna_start::Union{Int, Missing}
    rna_end::Union{Int, Missing}

    # Experimental data:
    expression::Union{Vector{Float64}, Missing}
end

struct Annotation <: Feature

    elements::Union{Dict{String, Vector{Feature}}, Missing} # Parent Node(s)
    annotation::String
end

struct Region <: Feature

    scaffold::Union{Scaffold, Missing} # Parent Node
    contig::Union{Contig, Missing} # Parent Node
    region_start::Union{Int, Missing}
    region_end::Union{Int, Missing}
    annotations::Union{Dict{String, Vector{Annotation}}, Missing}

    # Experimental data:
    signals::Union{Vector{Union{Vector{UInt8},Vector{UInt16},Vector{UInt32},Vector{UInt64}, Vector{Float64}}}, Nothing}
    binsignals::Union{Vector{BitVector}, Nothing}
    samples::Union{Vector{String}, Nothing}
end

struct Repeat <: Feature 

    scaffold::Union{Scaffold, Missing} # Parent Node
    contig::Union{Contig, Missing} # Parent Node
    regions::Union{Vector{Region}, Missing} 
    family::Union{String, Missing}
    type::Union{String, Missing}
    repeat_start::Union{Int, Missing}
    repeat_end::Union{Int, Missing}
    sequence::Union{B, Missing} where B <: BioSequences.BioSequence

    # Experimental data:
    signals::Union{Vector{Union{Vector{UInt8},Vector{UInt16},Vector{UInt32},Vector{UInt64}, Vector{Float64}}}, Nothing}
    binsignals::Union{Vector{BitVector}, Nothing}
    samples::Union{Vector{String}, Nothing}
end

struct Segment <: Feature

    prev_segment::Union{Vector{Feature}, Missing} 
    # plus_elements::Union{Dict{String, Vector{Feature}}, Missing} # Parent Node(s)
    # minus_elements::Union{Dict{String, Vector{Feature}}, Missing} # Parent Node(s)
    segment_start::Union{Int, Missing}
    segment_end::Union{Int, Missing}
    annotations::Union{Dict{String, Vector{Annotation}}, Missing}
    next_segment::Union{Vector{Feature}, Missing}

    # Experimental data:
    signals::Union{Vector{Union{Vector{UInt8},Vector{UInt16},Vector{UInt32},Vector{UInt64}, Vector{Float64}}}, Nothing}
    binsignals::Union{Vector{BitVector}, Nothing}
    samples::Union{Vector{String}, Nothing}
end

mutable struct RefGenome
    
    scaffolds::Dict{String, Scaffold}
    contigs::Dict{String, Contig}
    enhancers::Vector{Enhancer}
    promoters::Vector{Promoter}
    genes::Tuple{Vector{String}, Vector{Gene}}
    introns::Vector{Intron}
    exons::Vector{Exon}
    rnas::Vector{RNA}
    segments::Dict{String, Vector{Segment}}
    repeats::Vector{Repeat}
    regions::Vector{Region}
    annotations::Vector{Annotation}
end

abstract type RangeAnchor end
struct TSS <: RangeAnchor end
struct REGION <: RangeAnchor end
struct TES <: RangeAnchor end

struct GeneRange

    range_start::RangeAnchor
    range_stop::RangeAnchor
    start_offset::Int
    stop_offset::Int
end

# Constructors:
GeneRange(start::RangeAnchor, stop::RangeAnchor) = begin

    return GeneRange(start, stop, 0, 0)
end

# Empty constructor:
function RefGenome()

    return RefGenome(Dict{String, Scaffold}(), 
                      Dict{String, Contig}(), 
                      Vector{Enhancer}(), 
                      Vector{Promoter}(), 
                      (Vector{String}(), Vector{Gene}()), 
                      Vector{Intron}(), 
                      Vector{Exon}(), 
                      Vector{RNA}(), 
                      Dict{String, Vector{Segment}}(),
                      Vector{Repeat}(),
                      Vector{Region}(), 
                      Vector{Annotation}())
end

# Displaying custom types:

Base.show(io::IO, x::Feature) = begin

    type = typeof(x)
    fields = fieldnames(type)
    fields_missing = [ismissing(getfield(x, i)) for i in fields]
    fields_present = fields[.!(fields_missing)]
    fields_missing = fields[fields_missing]

    if any(:id .== fields_present)

        id = x.id

        print(io, "$type: '$id', with fields: ")
        for field in fields_present
            print(io, "'$field' ")
        end
    else

        print(io, "$type with fields: ")
        for field in fields_present
            print(io, "'$field' ")
        end
    end

    print(io, "\n     missing fields: ")
    for field in fields_missing
        print(io, "'$field' ")
    end
end

Base.show(io::IO, x::Vector{Feature}) = begin
    
    type = typeof(x)
    len = length(x)
    println(io, "$type with $len items")
end

export Scaffold,
       Contig,
       Junk,
       Enhancer,
       Gene,
       mRNA,
       Intron,
       Exon,
       Segment,
       Repeat,
       Annotation,
       Region,
       RefGenome