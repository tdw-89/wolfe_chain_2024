using Pkg
Pkg.activate(".")

# Need to add first, as it's unregistered
Pkg.add(url="https://github.com/tk3369/YeoJohnsonTrans.jl.git")

# Need to add my version of GFF3 for compat with automa v1
Pkg.add(url="https://github.com/tdw-89/GFF3.jl.git")

# Need to add my version of BED for narrowPeak compat
Pkg.add(url="https://github.com/tdw-89/BED.jl.git")

# Custom library
Pkg.add(url="https://github.com/tdw-89/BioinfoTools.git")

# List of external dependencies
external_deps = [
    "CSV", "DataFrames", "FastaIO", "Intervals", "IntervalTrees", "BioGenerics", "BioSequences", "BioAlignments",
    "Serialization", "StatsBase", "CategoricalArrays", "PlotlyJS", "Printf", "Pipe",
    "NaturalSort", "Interpolations", "XAM", "GenomicFeatures", "Indexes",
    "TranscodingStreams", "DataStructures", "Graphs", "MetaGraphs", "EzXML",
    "RollingFunctions", "DSP", "HypothesisTests", "MultipleTesting", "FASTX",
    "FreqTables", "Loess", "CodecZlib", "BGZFStreams",
    "FixedPointNumbers", "GLM", "Automa", "JSON", "Distributions", "StatsModels",
    "OhMyThreads", "GaussianMixtures", "LanguageServer"
]
Pkg.add(external_deps)
println("All packages installed")