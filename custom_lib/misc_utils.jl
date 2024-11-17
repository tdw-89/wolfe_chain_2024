#= An assortment of utility functions and corresponding types that don't fit anywhere else atm =#

if !@isdefined Feature
    include("genome_types.jl")

end

import Base.get

"""Check if two genes overlap"""
function hasoverlap(first_start::Int, second_start::Int, first_end::Int, second_end::Int)

    if (first_end <= second_start) || (second_end <= first_start)

        return false
    end

    return true
end

function hasoverlap(geneA::Gene, geneB::Gene)
    
    return hasoverlap(geneA.gene_start, geneB.gene_start, geneA.gene_end, geneB.gene_end)
end

"""Given two genes that overlap, return the overlap length"""
function overlaplength(first_start::Int, second_start::Int, first_end::Int, second_end::Int)

    if !hasoverlap(first_start, second_start, first_end, second_end)

        return 0
    end

    return max((min(first_end, second_end) - max(first_start, second_start)) + 1, 0)
end

function overlaplength(geneA::Gene, geneB::Gene)
    
    return overlaplength(geneA.gene_start, geneB.gene_start, geneA.gene_end, geneB.gene_end)
end

"""
    This function assumes that the 'paralog' DataFrame has (at least) the following columns in the
    following order, with the following types:
        * 1: GeneID <: AbstractString
        * 2: ParalogID <: AbstractString
    This function also assumes that the 'expression' DataFrame has (at least) the following columns in the
    following order, with the following types:
        * 1: GeneID <: AbstractString
        * 2: Expression value <: AbstractFloat
"""
function add_expression_to_paralogs(paralog_df::DataFrame, expression_df::DataFrame)

    gene_expr_vals = Vector{Union{Float64, Missing}}(undef, size(paralog_df)[1])
    paralog_expr_vals = Vector{Union{Float64, Missing}}(undef, size(paralog_df)[1])

    for i in 1:size(paralog_df)[1]

        gene_id = paralog_df[i,1]
        paralog_id = paralog_df[i,2]
        gene_expr_ind = findfirst(expression_df[:,1] .== gene_id)
        paralog_expr_ind = findfirst(expression_df[:,1] .== paralog_id)
        gene_expr_val = isnothing(gene_expr_ind) ? missing : expression_df[gene_expr_ind,2]
        paralog_expr_val = isnothing(paralog_expr_ind) ? missing : expression_df[paralog_expr_ind,2]

        gene_expr_vals[i] = gene_expr_val
        paralog_expr_vals[i] = paralog_expr_val
    end

    return hcat(paralog_df, DataFrame("GeneExpr" => gene_expr_vals, "ParalogExpr" => paralog_expr_vals))
end

function to_vector(gene_range::GeneRange)
    start_anchor = gene_range.range_start
    end_anchor = gene_range.range_stop

    if start_anchor != end_anchor

        error("Must supply a specific gene if start and end anchors for the GeneRange are different.")
    end

    return gene_range.start_offset:gene_range.stop_offset
end

function to_vector(gene_range::GeneRange, gene::Gene)
    error("not implemented yet!")

end

function sortgenes!(ref_genome::RefGenome)

    if isempty(ref_genome.scaffolds)

        @warn "No scaffold/chromosome information.\nTreating all genes as existing on a single chromosome"
        # TODO: implement
    else

        for scaffold in values(ref_genome.scaffolds)

            sort!(scaffold.genes, by=gene -> gene.gene_start)
        end
    end
end

# Overload 'get' function for RefGenome
function Base.get(genome::RefGenome, gene_id::S) where S <: AbstractString

    gene_ind = findfirst(genome.genes[1] .== gene_id)

    if isnothing(gene_ind)

        return missing
    else
            
        return genome.genes[2][gene_ind]
    end
end

function Base.get(genome::RefGenome, gene_ids::Vector{S}) where S <: AbstractString

    return_list = Vector{Union{Gene, Missing}}(undef, length(gene_ids))

    for (i,id) in enumerate(gene_ids)

        return_list[i] = get(genome, id)
    end

    return return_list
end 
