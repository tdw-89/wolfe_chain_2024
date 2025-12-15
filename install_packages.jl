using Pkg

full_list = [
    "CSV", "DataFrames", "FastaIO", "Intervals", "IntervalTrees", "BioGenerics", "BioSequences", "BioAlignments",
    "Serialization", "StatsBase", "CategoricalArrays", "PlotlyJS", "Printf", "Pipe",
    "NaturalSort", "Interpolations", "XAM", "GenomicFeatures", "Indexes",
    "TranscodingStreams", "DataStructures", "Graphs", "MetaGraphs", "EzXML",
    "RollingFunctions", "DSP", "HypothesisTests", "MultipleTesting", "FASTX",
    "GFF3", "FreqTables", "Loess", "CodecZlib", "BGZFStreams",
    "FixedPointNumbers", "GLM", "Automa"
]

to_install = String[]

for pkg in full_list

    try
        # Check if the package is already installed
        temp_expr = quote
            using $(Symbol(pkg))
        end
        eval(temp_expr)
    catch e
        push!(to_install, pkg)
    end
end

# Install missing packages
if !isempty(to_install)
    Pkg.add(to_install)
    
end

try
    temp_expr = quote
        using YeoJohnsonTrans
    end
    eval(temp_expr)
catch e
    Pkg.add("https://github.com/tk3369/YeoJohnsonTrans.jl")
end

println("All packages installed")
