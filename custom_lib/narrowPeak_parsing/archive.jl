mutable struct Archive
    num_records::Int
    records::Cons{Record}
end

"""
    PeakData.Archive()
Create an empty archive. A list with an empty record object will be placed
in the 'records' field.
"""
function Archive()
    return Archive(0,list(Record()))
end

function push!(archive::Archive,record::Record)
    if(archive.num_records == 0)
        archive.records = list(record)
        archive.num_records += 1
        nothing
    else
        archive.records = cons(record,archive.records)
        archive.num_records += 1
        nothing
    end
end

function Base.convert(archive::Archive) # Is overloading a base function OK?
    df = DataFrame(
    :chrom => fill("",archive.num_records),
    :start => fill(0,archive.num_records),
    :end => fill(0,archive.num_records),
    :peakName => fill("",archive.num_records),
    :score => fill(0,archive.num_records),
    :strand => fill('.',archive.num_records),
    :signalValue => fill(0.0,archive.num_records),
    :pValue => fill(0.0,archive.num_records),
    :qValue => fill(0.0,archive.num_records),
    :peakPoint => fill(0,archive.num_records),
    )
    i = archive.num_records
    for record in archive.records
        df.chrom[i] = chrom(record)
        df.start[i] = chromstart(record)
        df.end[i] = chromend(record)
        df.peakName[i] = name(record)
        df.score[i] = score(record)
        df.strand[i] = strand(record)
        df.signalValue[i] = signalValue(record)
        df.pValue[i] = pValue(record)
        df.qValue[i] = qValue(record)
        df.peakPoint[i] = peakPoint(record)
        i -= 1
    end

    return df
end