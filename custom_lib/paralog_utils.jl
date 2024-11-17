# TO DO:
# - Add a tie-breaking scheme to each scoring method
using DataFrames
using Graphs
using MetaGraphs
# using GraphIO
using EzXML

function rbh(paralog_df::DataFrame; scoring::String="max")

    scoring = lowercase(scoring)

    if !(scoring in ["max", "maximum", "double_max", "mean", "average", "avg"])

        error("Invalid scoring method. Must be 'max', 'maximum', 'double_max', 'mean', 'avg', or 'average'.")
    end

    scoring = scoring in ["max", "maximum"] ? "max" : scoring in ["mean", "avg", "average"] ? "mean" : scoring

    @assert typeof(paralog_df[1,1]) <: AbstractString
    @assert typeof(paralog_df[1,2]) <: AbstractString
    @assert typeof(paralog_df[1,3]) <: AbstractFloat
    @assert typeof(paralog_df[1,4]) <: AbstractFloat

    unique_ids = unique(vcat(paralog_df[:,1], paralog_df[:,2]))
    ids_to_ind_dict = Dict(unique_ids[i] => i for i in eachindex(unique_ids))
    ind_to_ids_dict = Dict(i => unique_ids[i] for i in eachindex(unique_ids))

    # Create a matrix of zeros
    rbh_matrix = zeros(Float64, length(unique_ids), length(unique_ids))
    orig_mat = zeros(Float64, length(unique_ids), length(unique_ids))

    # Fill in the matrix, treating gene 'i' as the query and gene 'j' as the subject
    for row in eachrow(paralog_df)

        i = ids_to_ind_dict[row[1]]
        j = ids_to_ind_dict[row[2]]
        
        orig_mat[i,j] = row[3]
        orig_mat[j,i] = row[4]

        if scoring == "max"

            max_perc_temp = max(row[3], row[4])
            rbh_matrix[i,j] = max_perc_temp
            rbh_matrix[j,i] = max_perc_temp
        elseif scoring == "mean"

            mean_perc_temp = (row[3] + row[4]) / 2
            rbh_matrix[i,j] = mean_perc_temp
            rbh_matrix[j,i] = mean_perc_temp
        else

            rbh_matrix[i,j] = orig_mat[i,j]
            rbh_matrix[j,i] = orig_mat[j,i]
        end
    end

    rbh_gene, rbh_paralog = String[], String[]
    matched_inds = Int[]
    perc_i_origs, perc_j_origs, max_percs, mean_percs = Float64[], Float64[], Float64[], Float64[]
        
    for i in 1:size(rbh_matrix)[1]

        perc_i, max_i = findmax(rbh_matrix[i,:])[1:2]
        perc_j, max_j = findmax(rbh_matrix[:,max_i])[1:2]
        perc_i_orig = orig_mat[i,max_i]
        perc_j_orig = orig_mat[max_i,i]
        max_perc = max(perc_i_orig, perc_j_orig)
        mean_perc = (perc_i_orig + perc_j_orig) / 2
        
        if i == max_j && max_i ∉ matched_inds && max_j ∉ matched_inds
            
            id_i = ind_to_ids_dict[max_i]
            id_j = ind_to_ids_dict[max_j]
            
            push!(rbh_gene, id_i)
            push!(rbh_paralog, id_j)
            push!(matched_inds, max_i)
            push!(matched_inds, max_j)
            push!(perc_i_origs, perc_i_orig)
            push!(perc_j_origs, perc_j_orig)
            push!(max_percs, max_perc)
            push!(mean_percs, mean_perc)
        end
    end

    return DataFrame("GeneID" => rbh_gene, "ParalogID" => rbh_paralog, "perc_1" => perc_i_origs, "perc_2" => perc_j_origs, "max_perc" => max_percs, "mean_perc" => mean_percs)
end

"""
    This function assumes that the input DataFrame has the following columns in the
    following order, with the following types:
        * 1: GeneID <: AbstractString
        * 2: ParalogID <: AbstractString
        * 3: edge value <: AbstractFloat (e.g. dS, dN, avg. % ID, etc.)
        * 4,5: (optionally) exprssion value for GeneID <: AbstractFloat and ParalogID <: AbstractFloat (respectively)
    This function also assumes that 'dS' is pairwise, so dS from a -> b is the same as dS from b -> a.
"""
function findfamilies(paralog_df::DataFrame; apply_cutoff::Bool=false,
                                            cutoff_val::Union{Float64, Nothing}, 
                                            cutoff_variable::Union{String, Symbol, Nothing}=nothing, 
                                            cutoff_comparison::F, 
                                            weighted_graph::Bool=false, 
                                            add_names::Bool=false) where F <: Function

    if ncol(paralog_df) != 5 && ncol(paralog_df) != 3

        error("Invalid input DataFrame. Must have 3, or 5 columns.")
    elseif ncol(paralog_df) == 5 && !(typeof(paralog_df[1,4]) <: AbstractFloat && typeof(paralog_df[1,5]) <: AbstractFloat)

        error("Invalid input DataFrame. Columns 4 and 5 must be <: AbstractFloat.")
    elseif typeof(paralog_df[1,1]) <: AbstractString && typeof(paralog_df[1,2]) <: AbstractString && typeof(paralog_df[1,3]) <: AbstractFloat

        nothing
    else

        error("Invalid input DataFrame. Items in columns 1 and 2 must be <: AbstractString, and column 3 must be <: AbstractFloat.")
    end

    if apply_cutoff
        
        filtered_df = filter(row -> cutoff_comparison(row[cutoff_variable], cutoff_val), paralog_df)
    end

    if weighted_graph

        return create_weighted_graph(filtered_df, cutoff_variable, cutoff_val=isnothing(cutoff_val) ? Inf : cutoff_val, add_names=add_names, add_expr=ncol(paralog_df) == 5) # Change this so adding exprsesion is an explicit parameter, not assumed (also change so that exprssion values can be in other columns)
    else

        return create_unweighted_graph(filtered_df, add_names=add_names)
    end
end

function family_sizes(name_vecs::Vector{Vector{String}})

    family_sizes = length.(name_vecs)
    gene_names, family_size = String[], Int[]
    
    for (i, vec) in enumerate(name_vecs)

        for gene_name in vec

            push!(gene_names, gene_name)
            push!(family_size, family_sizes[i])
        end
    end

    return DataFrame("GeneID" => gene_names, "FamilySize" => family_size)
end

function create_unweighted_graph(paralog_df::DataFrame; add_names::Bool=false)

    ids = String[]
    add_vert_1, add_vert_2 = true, true
    dup_graph = SimpleGraph()
    
    for i in 1:size(paralog_df)[1]

        id_1 = paralog_df[i,1]
        id_2 = paralog_df[i,2]

        id_1_vert = findfirst(ids .== id_1)
        id_2_vert = findfirst(ids .== id_2)

        if isnothing(id_1_vert)

            push!(ids, id_1)
            id_1_vert = length(ids)
            add_vert_1 = true
        else

            add_vert_1 = false
        end

        if isnothing(id_2_vert)

            push!(ids, id_2)
            id_2_vert = length(ids) 
            add_vert_2 = true
        else

            add_vert_2 = false
        end

        if add_vert_1
                
            add_vertex!(dup_graph)
        end

        if add_vert_2

            add_vertex!(dup_graph)
        end

        if has_edge(dup_graph, id_1_vert, id_2_vert)

            nothing
        else

            add_edge!(dup_graph, id_1_vert, id_2_vert)
        end
    end

    if add_names

        dup_graph = MetaGraph(dup_graph, 1)
        
        for i in eachindex(ids)

            set_prop!(dup_graph, i, :name, ids[i])
        end
    end

    return dup_graph
end

function create_weighted_graph(paralog_df::DataFrame, weight_variable::Union{String, Symbol}; cutoff_val::Float64=Inf, add_names::Bool=false, add_expr::Bool=false)

    ids = String[]
    add_vert_1, add_vert_2 = true, true
    dup_graph = MetaGraph(SimpleGraph(), cutoff_val)

    for i in 1:size(paralog_df)[1]

        id_1 = paralog_df[i,1]
        id_2 = paralog_df[i,2]

        id_1_vert = findfirst(ids .== id_1)
        id_2_vert = findfirst(ids .== id_2)

        if isnothing(id_1_vert)

            push!(ids, id_1)
            id_1_vert = length(ids)
            add_vert_1 = true
        else

            add_vert_1 = false
        end

        if isnothing(id_2_vert)

            push!(ids, id_2)
            id_2_vert = length(ids) 
            add_vert_2 = true
        else

            add_vert_2 = false
        end

        if add_vert_1

            add_vertex!(dup_graph)
        end

        if add_vert_2

            add_vertex!(dup_graph)
        end

        if has_edge(dup_graph, id_1_vert, id_2_vert)
            
            nothing
        else

            add_edge!(dup_graph, id_1_vert, id_2_vert)
            set_prop!(dup_graph, id_1_vert, id_2_vert, :weight, paralog_df[i, weight_variable])

            if add_expr
                    
                set_prop!(dup_graph, id_1_vert, :expr, paralog_df[i,4])
                set_prop!(dup_graph, id_2_vert, :expr, paralog_df[i,5])
            end
        end
    end

    if add_names

        for i in 1:length(ids)

            set_prop!(dup_graph, i, :name, ids[i])
        end
    end

    return dup_graph
end

function savegraphcsv(graph::MG, filename::String; force::Bool=false, add_names::Bool=false) where MG <: MetaGraph

    if filename in readdir() && !force

        error("File already exists. Use 'force=true' to overwrite.")
    end

    open(filename, "w") do io

        if add_names

            println(io, "Source,Target,Weight,Source_Name,Target_Name")
        else
           
            println(io, "Source,Target,Weight")
        end

        for (edge,weight_dict) in graph.eprops

            temp_str_vec = split(string(edge), " => ")
            v1 = match(r"[0-9]+", temp_str_vec[1]).match
            v2 = match(r"[0-9]+", temp_str_vec[2]).match

            if add_names
                
                v1_name = props(graph, parse(Int, v1))[:name]
                v2_name = props(graph, parse(Int, v2))[:name]

                println(io, "$v1_name,$v2_name,$(10 ^ weight_dict[:weight]),$v1_name,$v2_name")
            else
                
                println(io, "$v1,$v2,$(10 ^ weight_dict[:weight])")
            end
        end
    end
end
