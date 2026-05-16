using DataFrames, CSV, OhMyThreads

const λ = 1/20_000

function calc_density_scores(line)
    cols = split(strip(line), '\t')

    gene_id = cols[4]
    strand  = cols[6]

    gene_start  = parse(Int, cols[2])
    gene_end    = parse(Int, cols[3])
    gene_length = gene_end - gene_start

    te_start = parse(Int, cols[8])
    te_end   = parse(Int, cols[9])

    int_start = max(te_start, gene_start)
    int_end   = min(te_end, gene_end)
    internal_score = max(0, int_end - int_start) / gene_length

    left_score = 0.0
    for p in te_start:(min(te_end, gene_start) - 1)
        distance = gene_start - p
        left_score += 1.0 / exp(distance * λ)
    end

    right_score = 0.0
    for p in max(te_start, gene_end):(te_end - 1)
        distance = p - gene_end
        right_score += 1.0 / exp(distance * λ)
    end

    upstream_score, downstream_score = strand == "+" ? (left_score, right_score) : (right_score, left_score)

    return gene_id, upstream_score, internal_score, downstream_score
end

function build_dict(id_score_tuples)
    dict = Dict{String, NTuple{3, Float64}}()
    for (gene_id, upstream, internal, downstream) in id_score_tuples
        if haskey(dict, gene_id)
            u, i, d = dict[gene_id]
            dict[gene_id] = (u + upstream, i + internal, d + downstream)
        else
            dict[gene_id] = (upstream, internal, downstream)
        end
    end
    return dict
end

function create_df(gene_scores::Dict{String, NTuple{3, Float64}})
    gene_ids   = collect(keys(gene_scores))
    upstreams  = [gene_scores[g][1] for g in gene_ids]
    internals  = [gene_scores[g][2] for g in gene_ids]
    downstreams = [gene_scores[g][3] for g in gene_ids]
    return DataFrame(
        Gene_ID = gene_ids,
        Upstream = upstreams, 
        Internal = internals, 
        Downstream = downstreams
        )
end

human = true

window_file = human ? "../../dicty_data/te_density/tes_in_window_human.bed" : "../../dicty_data/te_density/tes_in_window.bed"

gene_scores = Dict{String, Float64}()
fh = open(window_file, "r")
lines = readlines(fh)
close(fh)

df = map(calc_density_scores, lines) |> build_dict |> create_df

output_file = human ? "../../dicty_data/te_density/te_density_human.csv" : "../../dicty_data/te_density/te_density.csv"
CSV.write(output_file, df)