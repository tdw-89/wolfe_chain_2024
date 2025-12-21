include("prelude.jl")

re_load_interact_data = false

# STRING protein-protein interaction data:
prot_interact_file = "../../dicty_data/mammals/primates/h_sapiens/STRING/9606.protein.links.full.v12.0.txt"
prot_alias_file = "../../dicty_data/mammals/primates/h_sapiens/STRING/9606.protein.aliases.v12.0.txt"

# Load the protein-protein interaction data:
prot_interact_df = CSV.read(prot_interact_file, DataFrame)

prot_alias_df = CSV.read(prot_alias_file, DataFrame)

if re_load_interact_data

    # Filter the protein-protein interaction data to high-confidence interactions:
    # NOTE: Right now we are just using the combined score, but as a lot of interactions
    # seem to come from text-mining, this data source type should be researched more.
    prot_interact_df = filter(row -> row.combined_score >= 700, prot_interact_df)

    # Filter to unique entries:
    sets_col = [Set([row.protein1, row.protein2]) for row in eachrow(prot_interact_df)]
    prot_interact_df = hcat(DataFrame(:Set => sets_col), select(prot_interact_df, Not([1,2])))
    unique!(prot_interact_df)
    set_list = collect.(prot_interact_df.Set)
    prot_1 = [set[1] for set in set_list]
    prot_2 = [set[2] for set in set_list]
    prot_interact_df = hcat(DataFrame(:protein1 => prot_1, :protein2 => prot_2), select(prot_interact_df, Not(1)))


    # Rename the STRING protein IDs to gene IDs:
    filter!(row -> contains(row[2], "ENSG"), prot_alias_df)
    prot_alias_df = unique(select(prot_alias_df, 1, 2))
    alias_dict = Dict(pair[1] => pair[2] for pair in zip(prot_alias_df[!,1], prot_alias_df[!,2]))
    prot_interact_df.protein1 = map(id -> get(alias_dict, id, id), prot_interact_df.protein1)
    prot_interact_df.protein2 = map(id -> get(alias_dict, id, id), prot_interact_df.protein2)

    # Save the protein-protein interaction data:
    CSV.write("../../dicty_data/human_prot_interact.csv", prot_interact_df)
else

    prot_interact_df = CSV.read("../../dicty_data/human_prot_interact.csv", DataFrame)
end

# Create a protein interaction count dictionary:
unique_genes = unique(vcat(prot_interact_df.protein1, prot_interact_df.protein2))
prot_interact_count = Dict(gene => 0 for gene in unique_genes)
for gene in vcat(prot_interact_df.protein1, prot_interact_df.protein2)
    prot_interact_count[gene] += 1
end

# Save the protein interaction count dictionary:
count_df = DataFrame(:gene => collect(keys(prot_interact_count)), :count => collect(values(prot_interact_count)))
CSV.write("../../dicty_data/human_prot_interact_count.csv", count_df)