module PeakData

import Automa
import Automa.RegExp: @re_str
import BGZFStreams
import BioGenerics
import FixedPointNumbers: N0f8

using GenomicFeatures
using Indexes
using TranscodingStreams
using DataFrames
using DataStructures

include("record.jl")
include("archive.jl")
include("reader.jl")
include("writer.jl")

function LoadNarrowPeak(path::AbstractString)
    if !occursin(".narrowPeak", path) && !occursin(".bed", path)
        e = "A .narrowPeak file was not specified by <path> argument"
        throw(e)
    end
    
    reader = open(PeakData.Reader, path)
    archive = PeakData.Archive()
    for record in reader
        PeakData.push!(archive, record)
    end
    close(reader)
    return convert(archive)
end

export LoadNarrowPeak

end # module

