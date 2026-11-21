using DataFrames, CSV

# Distance decay: 4.5 kb, the mean extent of H3K9me2 spread from euchromatic TEs in
# D. melanogaster (Lee & Karpen 2017, eLife 6:e25762; mean 4.50 kb, SD 4.59 kb, 1 kb windows).
const λ = 1/4_500

# Flank window used to build tes_in_window*.bed (bedtools window -w 100000).
const WINDOW = 100_000

# Score of a flank that is entirely TE: sum of exp(-λd) for d = 1..WINDOW. Dividing by it puts
# Upstream/Downstream on a 0-1 scale (distance-weighted fraction of the flank covered by TEs), the same
# scale as Internal (fraction of the gene covered by TEs), so the three can be compared or summed.
# Genes closer than WINDOW to a chromosome end are not corrected for their shorter flank.
const FLANK_MAX = sum(exp(-λ * d) for d in 1:WINDOW)

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

    # TE bases are clipped to the window so a fully TE-covered flank scores exactly 1
    left_score = 0.0
    for p in max(te_start, gene_start - WINDOW):(min(te_end, gene_start) - 1)
        distance = gene_start - p        # 1 for the base next to the gene
        left_score += exp(-distance * λ)
    end

    right_score = 0.0
    for p in max(te_start, gene_end):(min(te_end, gene_end + WINDOW) - 1)
        distance = p - gene_end + 1      # 1 for the base next to the gene
        right_score += exp(-distance * λ)
    end

    left_score  /= FLANK_MAX
    right_score /= FLANK_MAX

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

# Both species; the window files are written by te_windows.sh
te_density_dir = joinpath(@__DIR__, "../../dicty_data/te_density")
for tag in ["", "_human"]
    lines = readlines(joinpath(te_density_dir, "tes_in_window$tag.bed"))
    df = map(calc_density_scores, lines) |> build_dict |> create_df
    CSV.write(joinpath(te_density_dir, "te_density$tag.csv"), df)
end
