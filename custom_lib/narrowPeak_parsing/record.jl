# Adapted from BED.jl source code (Ciarán O'Mara), specifically
# the record.jl file.

# NOTES:
#       > Sometimes the original program checked to see if a record contained a specified 
#           value just by checking the number of columns. This may lead to inflexibility.
#           For an example see 'hasqValue' function. Same thing for functions that invoke 'isfilled()'.
#       > Check the (ensembl) rules for the peak point-source ('peakPoint') before creating
#           the 'peakPoint' function. you noted that the point-source was 0-based in the
#           doc-string for its accessor, is this correct?


mutable struct Record
    # data and filled range
    data::Vector{UInt8}
    filled::UnitRange{Int}
    # number of columns
    ncols::Int
    # indexes
    chrom::UnitRange{Int}
    chromstart::UnitRange{Int}
    chromend::UnitRange{Int}
    name::UnitRange{Int}
    score::UnitRange{Int}
    strand::Int
    signalValue::UnitRange{Int} 
    pValue::UnitRange{Int} 
    qValue::UnitRange{Int} 
    peakPoint::UnitRange{Int} 
end

"""
    PeakData.Record()
Create an unfilled peak data record.
"""
function Record()
    return Record(
        UInt8[], 1:0, 0,
        # chrom-score
        1:0, 1:0, 1:0, 1:0, 1:0,
        # strand-qValue
        0, 1:0, 1:0, 1:0,
        # peakPoint
        1:0)
end

"""
PeakData.Record(data::Vector{UInt8})
Create a eak data record object from `data`.
This function verifies and indexes fields for accessors.
Note that the ownership of `data` is transferred to a new record object.
"""
function Record(data::Vector{UInt8})
    return convert(Record, data)
end

function Base.convert(::Type{Record}, data::Vector{UInt8})
    record = Record(
        data, 1:0, 0,
        # chrom-score
        1:0, 1:0, 1:0, 1:0, 1:0,
        # strand-qValue
        0, 1:0, 1:0, 1:0,
        # peakPoint
        1:0)
    index!(record)
    return record
end

"""
    BED.Record(str::AbstractString)
Create a BED record object from `str`.
This function verifies and indexes fields for accessors.
"""
function Record(str::AbstractString)
    return convert(Record, str)
end

function Base.convert(::Type{Record}, str::AbstractString)
    return convert(Record, Vector{UInt8}(str))
end

function Base.empty!(record::Record)
    record.filled = 1:0
    record.ncols = 0
    record.chrom = 1:0
    record.chromstart = 1:0
    record.chromend = 1:0
    record.name = 1:0
    record.score = 1:0
    record.strand = 0
    record.signalValue = 1:0
    record.pValue = 1:0
    record.qValue = 1:0
    record.peakPoint = 1:0
    return record
end

# This function returns an 'Interval' object created from the 'Record' object.
# The functions from the 'BioGenerics' library are overloaded further down this
# file in the 'Accessor' functions section. 'strand' and 'hasstrand' are also 
# defined further down.
function GenomicFeatures.Interval(record::Record)
    name = BioGenerics.seqname(record)
    lpos = BioGenerics.leftposition(record)
    rpos = BioGenerics.rightposition(record)
    strd = hasstrand(record) ? GenomicFeatures.strand(record) : GenomicFeatures.STRAND_BOTH
    return GenomicFeatures.Interval(name, lpos, rpos, strd, record)
end

function Base.convert(::Type{GenomicFeatures.Interval}, record::Record)+ 1
    return GenomicFeatures.Interval(record)
end

function Base.convert(::Type{GenomicFeatures.Interval{Record}}, record::Record)
    return convert(GenomicFeatures.Interval, record)
end

function isfilled(record::Record)
    return !isempty(record.filled)
end

function Base.:(==)(record1::Record, record2::Record)
    if isfilled(record1) == isfilled(record2) == true
        r1 = record1.filled
        r2 = record2.filled
        return length(r1) == length(r2) && memcmp(pointer(record1.data, first(r1)), pointer(record2.data, first(r2)), length(r1)) == 0
    end

    return isfilled(record1) == isfilled(record2) == false
end

function Base.copy(record::Record)
    return Record(
        record.data[record.filled],
        record.filled,
        record.ncols,
        record.chrom,
        record.chromstart,
        record.chromend,
        record.name,
        record.score,
        record.strand,
        record.signalValue,
        record.pValue,
        record.qValue,
        record.peakPoint)
end

function Base.write(io::IO, record::Record)
    return unsafe_write(io, pointer(record.data, first(record.filled)), length(record.filled)) # NEED TO CHANGE ?, NOT APPLICABLE FUNCTION
end

function Base.print(io::IO, record::Record)
    write(io, record)
    return nothing
end

function Base.show(io::IO, record::Record)
    print(io, summary(record), ':')
    if isfilled(record)
        println(io)
        println(io, "   chromosome: ", chrom(record))
        println(io, "        start: ", chromstart(record))
        print(io, "          end: ", chromend(record))
        if hasname(record)
            println(io)
            print(io, "         name: ", name(record))
        end
        if hasscore(record)
            println(io)
            print(io, "        score: ", score(record))
        end
        if hasstrand(record)
            println(io)
            print(io, "       strand: ", strand(record))
        end
        if hassignalValue(record)
            println(io)
            print(io, "  enrichment (signal value): ", signalValue(record))
        end
        if haspValue(record)
            println(io)
            print(io, "    p-value: ", pValue(record))
        end
        if hasqValue(record)
            println(io)
            print(io, "     q-value: ", qValue(record))
        end
        if haspeakPoint(record)
            println(io)
            print(io, "  peak point-source: ", peakPoint(record))
        end
    else
        print(io, " <not filled>")
    end
end


# Accessor functions
# ------------------

"""
    chrom(record::Record)::String
Get the chromosome name of `record`.
"""
function chrom(record::Record)::String
    checkfilled(record)
    return String(record.data[record.chrom])
end

function haschrom(record::Record)
    return isfilled(record)
end

function BioGenerics.seqname(record::Record)
    return chrom(record)
end

function BioGenerics.hasseqname(record::Record)
    return haschrom(record)
end

"""
    chromstart(record::Record)::Int
Get the starting position of `record`.
Note that the first base is numbered 1.
"""
function chromstart(record::Record)::Int
    checkfilled(record)
    return unsafe_parse_decimal(Int, record.data, record.chromstart) # MAYBE NEED TO CHANGE? Removed the + 1, as .narrowPeak files start on 0? 
end

function haschromstart(record::Record)
    return isfilled(record)
end

function BioGenerics.leftposition(record::Record)
    return chromstart(record)
end

function BioGenerics.hasleftposition(record::Record)
    return haschromstart(record)
end

"""
    chromend(record::Record)::Int
Get the end position of `record`.
"""
function chromend(record::Record)::Int
    checkfilled(record)
    return unsafe_parse_decimal(Int, record.data, record.chromend) # NEED TO CHANGE?, NOT APPLICABLE FUNCTION
end

function haschromend(record::Record)
    return isfilled(record)
end

function BioGenerics.rightposition(record::Record)
    return chromend(record)
end

function BioGenerics.hasrightposition(record::Record)
    return haschromend(record)
end

"""
    name(record::Record)::String
Get the name of `record`.
"""
function name(record::Record)::String
    checkfilled(record)
    if !hasname(record)
        missingerror(:name)
    end
    return String(record.data[record.name])
end

function hasname(record::Record)
    return record.ncols ≥ 4
end

"""
    score(record::Record)::Int
Get the score between 0 and 1000.
"""
function score(record::Record)::Int
    checkfilled(record)
    if !hasscore(record)
        missingerror(:score)
    end
    return unsafe_parse_decimal(Int, record.data, record.score) 
end

function hasscore(record::Record)
    return record.ncols ≥ 5
end

"""
    strand(record::Record)::GenomicFeatures.Strand
Get the strand of `record`.
"""
function strand(record::Record)::GenomicFeatures.Strand
    checkfilled(record)
    if !hasstrand(record)
        missingerror(:strand)
    end
    return convert(GenomicFeatures.Strand, Char(record.data[record.strand]))
end

function hasstrand(record::Record)
    return record.ncols ≥ 6
end

function GenomicFeatures.strand(record::Record)
    return strand(record)
end

"""
    signalValue(record::Record)::Float64
Get the enrichment value for the peak specified in 'record'.
"""
function signalValue(record::Record)::Float64
    checkfilled(record)
    if !hassignalValue(record)
        missingerror(:signalValue)
    end
    return parse_float(record.data, record.signalValue)# NEED TO CHANGE, Temporary solution (slow i think)
end

function hassignalValue(record::Record)
    return record.ncols ≥ 7
end

"""
    pValue(record::Record)::Float64
Get the p-value for the peak specified in 'record'.
"""
function pValue(record::Record)::Float64
    checkfilled(record)
    if !haspValue(record)
        missingerror(:pValue)
    end
    return parse_float(record.data, record.pValue) # NEED TO CHANGE, Temporary solution (slow i think)
end

function haspValue(record::Record)
    return record.ncols ≥ 8
end

"""
    qValue(record::Record)::Float64
Get the q-value for the peak specified in 'record'.
"""
function qValue(record::Record)::Float64
    checkfilled(record)
    if !hasqValue(record)
        missingerror(:qValue)
    end
    return parse_float(record.data, record.qValue) # NEED TO CHANGE, Temporary solution (slow i think)
end

function hasqValue(record::Record)
    return record.ncols ≥ 9
end

# function parse_rgbcolor(data::Vector{UInt8}, range::UnitRange{Int}) # DELETE ?
#     function searchcomma(s)
#         i = s
#         while i ≤ last(range)
#             if data[i] == UInt8(',')
#                 break
#             end
#             i += 1
#         end
#         return i
#     end
#     lo = first(range)
#     hi = searchcomma(lo)
#     r = unsafe_parse_byte(data, lo:hi - 1)
#     if hi > last(range)
#         # single value
#         g = b = r
#     else
#         # triplet
#         lo = hi + 1
#         hi = searchcomma(lo)
#         g = unsafe_parse_byte(data, lo:hi - 1)
#         lo = hi + 1
#         hi = searchcomma(lo)
#         b = unsafe_parse_byte(data, lo:hi - 1)
#     end
#     return ColorTypes.RGB(reinterpret(N0f8, r), reinterpret(N0f8, g), reinterpret(N0f8, b))
# end

# function unsafe_parse_byte(data::Vector{UInt8}, range::UnitRange{Int}) # DELETE ?
#     val::UInt8 = 0x00
#     for i in range
#         val = val * 0x0a + (data[i] - UInt8('0'))
#     end
#     return val
# end

"""
    peakPoint(record::Record)::UnitRange{Int}
Get the point-source for the peak specified in 'record' (1-based).
"""
function peakPoint(record::Record)::Int
    checkfilled(record)
    if !haspeakPoint(record)
        missingerror(:peakPoint)
    end
    return unsafe_parse_decimal(Int, record.data, record.peakPoint) # MAYBE NEED TO CHANGE ?
end

function haspeakPoint(record::Record)
    return record.ncols ≥ 10
end

function checkfilled(record::Record)
    if !isfilled(record)
        throw(ArgumentError("unfilled peak record"))
    end
end

# r"[-+]?[0-9]+" must match `data[range]`.
function unsafe_parse_decimal(::Type{T}, data::Vector{UInt8}, range::UnitRange{Int}) where T <: Signed
    lo = first(range)
    if data[lo] == UInt8('-')
        sign = T(-1)
        lo += 1
    elseif data[lo] == UInt8('+')
        sign = T(+1)
        lo += 1
    else
        sign = T(+1)
    end
    x = zero(T)
    @inbounds for i in lo:last(range)
        x = Base.Checked.checked_mul(x, 10 % T)
        x = Base.Checked.checked_add(x, (data[i] - UInt8('0')) % T)
    end
    return sign * x
end

function parse_float(data::Vector{UInt8}, range::UnitRange{Int}) # NEED TO CHANGE, Temporary solution (slow i think)
    float_string = Char.(data[range])
    float_string = join(float_string)
    return tryparse(Float64,float_string)
end

function memcmp(p1::Ptr, p2::Ptr, n::Integer)
    return ccall(:memcmp, Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t), p1, p2, n)
end
