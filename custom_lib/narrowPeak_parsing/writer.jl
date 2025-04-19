# Adapted from BED.jl source code (Ciarán O'Mara), specifically
# the writer.jl file.

"""
    BED.Writer(output::IO)
Create a data writer of the peak data.
# Arguments:
* `output`: data sink
"""
struct Writer <: BioGenerics.IO.AbstractWriter
    output::IO
end

function BioGenerics.IO.stream(writer::Writer)
    return writer.output
end

function Base.write(writer::Writer, record::Record)
    return write(writer.output, record, '\n')
end