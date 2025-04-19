# Adapted from BED.jl source code (Ciarán O'Mara), specifically
# the reader.jl file.
import Automa
import Automa.RegExp: @re_str
import Automa.Stream: @mark, @markpos, @relpos, @abspos

# Not sure what this does:
function appendfrom!(dst, dpos, src, spos, n)
    if length(dst) < dpos + n - 1
        resize!(dst, dpos + n - 1)
    end
    unsafe_copyto!(dst, dpos, src, spos, n)
    return dst
end

mutable struct Reader <: BioGenerics.IO.AbstractReader
    state::BioGenerics.Automa.State
    index::Union{Indexes.Tabix,Nothing}

    # Constructor:
    function Reader(stream::TranscodingStream, index=nothing)
        # I think this says to create an Automa state using 'stream'(arg 1)
        # with a starting state of 1 (arg 2), a starting line number
        # of 1, and unfilled status (arg 4). The call to 'new()' is
        # only used in the context of an object in a constructor method.
        return new(BioGenerics.Automa.State(stream, 1, 1, false), index)  
    end
end

# Function for creating a 'Reader' object given a 
# file stream (arg 1), and possibly an index file
# path (arg 2).
function Reader(input::IO; index=nothing)
    if isa(index, AbstractString)
        index = Indexes.Tabix(index)
    end

    stream = TranscodingStreams.NoopStream(input)

    # Call the 'Reader' object constructor. I think the dispatcher
    # sees the 'stream' argument and discerns that you want the constructor,
    # it is not a recursive call to this function itself.
    return Reader(stream, index)

end

# Version of 'Reader' object creating function that accepts a file path,
# instead of a file stream as input.
function Reader(filepath::AbstractString; index=:auto)
    if isa(index, Symbol) && index != :auto
        throw(ArgumentError("invalid index argument: ':$(index)'"))
    end
    if endswith(filepath, ".bgz")
        input = BGZFStreams.BGZFStream(filepath)
        if index == :auto
            index = Indexes.findtabix(filepath)
        end
    else
        input = open(filepath)
    end
    return Reader(input, index = index)
end

# Tell the function 'Base.eltype()' to return the 'Record'
# type when querying the element type of a 'Reader' iterator.
# In other words, if you create a for-loop using a 'Reader', then
# ask what type you are iterating through using 'eltype(Reader)',
# it will return INCOMPLETE
function Base.eltype(::Type{Reader})
    return Record
end

# Return the underlying file stream of the 'Reader' object
function BioGenerics.IO.stream(reader::Reader)
    return reader.state.stream
end

function Base.iterate(reader::Reader, nextone = Record())
    if BioGenerics.IO.tryread!(reader, nextone) === nothing
        return nothing
    end
    return copy(nextone), empty!(nextone) # Empty record for inplace reading and reuse of array allocations.
end

function GenomicFeatures.eachoverlap(reader::Reader, interval::GenomicFeatures.Interval)
    if reader.index === nothing
        throw(ArgumentError("index is null"))
    end
    return Indexes.TabixOverlapIterator(reader, interval)
end

const record_machine, file_machine = (function ()
    cat = Automa.RegExp.cat
    rep = Automa.RegExp.rep
    opt = Automa.RegExp.opt

    record = let
        chrom = re"[ -~]*" # NEED TO CHANGE? why isn't it [A-Za-z_]? why is a space ⎵ included as a possible character?
        chrom.actions[:enter] = [:pos]
        chrom.actions[:exit] = [:record_chrom]

        chromstart = re"[0-9]+"
        chromstart.actions[:enter] = [:pos]
        chromstart.actions[:exit] = [:record_chromstart]

        chromend = re"[0-9]+"
        chromend.actions[:enter] = [:pos]
        chromend.actions[:exit] = [:record_chromend]

        name = re"[ -~]*" # NEED TO CHANGE? why isn't it [A-Za-z_]? why is a space ⎵ included as a possible character?
        name.actions[:enter] = [:pos]
        name.actions[:exit] = [:record_name]

        score = re"[0-9]+"
        score.actions[:enter] = [:pos]
        score.actions[:exit] = [:record_score]

        strand = re"[+\-.?]"
        strand.actions[:enter] = [:record_strand] #Note: single byte.

        signalValue = re"[+\-0-9e\.]+" # NOT SURE if this is the correct RegExp
        signalValue.actions[:enter] = [:pos]
        signalValue.actions[:exit] = [:record_signalValue]

        pValue = re"[+\-0-9e\.]+" # NOT SURE if this is the correct RegExp
        pValue.actions[:enter] = [:pos]
        pValue.actions[:exit] = [:record_pValue]

        qValue = re"[+\-0-9e\.]+" # NOT SURE if this is the correct RegExp
        qValue.actions[:enter] = [:pos]
        qValue.actions[:exit] = [:record_qValue]

        peakPoint = re"[0-9]+"
        peakPoint.actions[:enter] = [:pos]
        peakPoint.actions[:exit] = [:record_peakPoint]

        cat(
            chrom, '\t',
            chromstart, '\t',
            chromend,
            opt(cat('\t', name,
            opt(cat('\t', score,
            opt(cat('\t', strand,
            opt(cat('\t', signalValue,
            opt(cat('\t', pValue,
            opt(cat('\t', qValue,
            opt(cat('\t', peakPoint)))))))))))))))
    end
    record.actions[:enter] = [:mark]
    record.actions[:exit] = [:record]

    newline = let
        lf = re"\n"
        lf.actions[:enter] = [:countline]

        cat(opt('\r'), lf)
    end

    file = rep(cat(record, newline))

    return map(Automa.compile, (record, file))
end)()


#= write("peak.dot", Automa.machine2dot(file_machine))
run(`dot -Tsvg -o peak.svg peak.dot`) =#

const record_actions = Dict(
    :mark => :(@mark),
    :pos => :(pos = @relpos(p)),
    :countline => :(),
    :record_chrom => :(record.chrom = (pos:@relpos(p-1)); record.ncols += 1),
    :record_chromstart => :(record.chromstart = (pos:@relpos(p-1)); record.ncols += 1),
    :record_chromend => :(record.chromend = (pos:@relpos(p-1)); record.ncols += 1),
    :record_name => :(record.name = (pos:@relpos(p-1)); record.ncols += 1),
    :record_score => :(record.score = (pos:@relpos(p-1)); record.ncols += 1),
    :record_strand => :(record.strand = @relpos(p); record.ncols += 1),
    :record_signalValue => :(record.signalValue = (pos:@relpos(p-1)); record.ncols += 1),
    :record_pValue => :(record.pValue = (pos:@relpos(p-1)); record.ncols += 1),
    :record_qValue => :(record.qValue = (pos:@relpos(p-1)); record.ncols += 1),
    :record_peakPoint => :(record.peakPoint = (pos:@relpos(p-1)); record.ncols += 1),
    :record => :(record.filled = 1:@relpos(p-1))
)

# Why are certain fields commented out and unfinished?
Automa.Stream.generate_reader(
    :index!,
    record_machine,
    arguments = (:(record::Record),),
    actions = record_actions,
    # context = :(),
    initcode = :(pos = 0),
    # loopcode = :()
    # returncode = :()
) |> eval


const initcode = quote
    pos = 0
    linenum = 0
    found_record=false
    # empty!(record)
    cs, linenum = state
end

const loopcode = quote
    if found_record
        @goto __return__
    end
end

Automa.Stream.generate_reader(
    :readrecord!,
    file_machine,
    arguments = (:(record::Record), :(state::Tuple{Int,Int})),
    actions = merge(record_actions, Dict(
        :record => quote
            appendfrom!(record.data, 1, data, @markpos, p-@markpos)
            record.filled = 1:(p-@markpos)
            found_record = true
            @escape
        end,
        :countline => :(linenum += 1),
    )),
    initcode = initcode,
    loopcode = loopcode,
    returncode = :(return cs, linenum, found_record)
) |> eval


function index!(record::Record)
    stream = TranscodingStreams.NoopStream(IOBuffer(record.data))
    cs = index!(stream, record)
    if cs < 0
        throw(ArgumentError("invalid record"))
    end
    return record
end

"""
    read!(rdr::Reader, rec::Record)
Read a `Record` into `rec`; overwriting or adding to existing field values.
It is assumed that `rec` is already initialized or empty.
"""
function Base.read!(rdr::Reader, record::Record)

    cs, ln, found = readrecord!(rdr.state.stream, record, (rdr.state.state, rdr.state.linenum))

    rdr.state.state = cs
    rdr.state.linenum = ln
    rdr.state.filled = found

    if found
        return record
    end

    if cs == 0 || eof(rdr.state.stream)
        throw(EOFError())
    end

    throw(ArgumentError("malformed file"))
end