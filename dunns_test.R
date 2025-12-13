lslibrary(dunn.test)

dunns_test <- function(means_vec) {
  means_vec |> dunn.test(method = "holm") -> d.test
  sig_comps <- d.test$P.adjusted < 0.05
  sig_comp_labels <- d.test$comparisons[sig_comps]
  sig_zs <- d.test$Z[sig_comps]
  sig_ps <- d.test$P.adjusted[sig_comps]
  return(
    list(
      "test.result" = d.test,
      "sig.comparisons" = sig_comp_labels,
      "sig.zs" = sig_zs,
      "sig.ps" = sig_ps
    )
  )
}

means_vecs_file <- "/Users/thomasd.wolfe/Library/CloudStorage/GoogleDrive-eccwolfe@gmail.com/My Drive/UML/Grad_school/research/code/analysis_scripts/dicty/final_code/data/julia_serialized/means_vecs_ds.json"
means_vecs_human_file <- "/Users/thomasd.wolfe/Library/CloudStorage/GoogleDrive-eccwolfe@gmail.com/My Drive/UML/Grad_school/research/code/analysis_scripts/dicty/final_code/data/julia_serialized/means_vecs_dS_human.json"
expr_deciles <- "/Users/thomasd.wolfe/Library/CloudStorage/GoogleDrive-eccwolfe@gmail.com/My Drive/UML/Grad_school/research/code/analysis_scripts/dicty/final_code/data/julia_serialized/expression_deciles.json"

# Enrich vs. dS:
vec_list <- jsonlite::fromJSON(means_vecs_file)
names(vec_list$K4me3) <- as.character(1:10)
vec_list$K4me3 |> dunns_test() -> k4.dunn.test
write.csv(as.data.frame(k4.dunn.test$test.result), "../../dicty_data/k4_dunns_test.csv")

# Expression vs. dS:
vec_list_expr <- jsonlite::fromJSON(expr_deciles)
dup_list <- vec_list_expr[-1]
dup_list <- dup_list[order(as.integer(names(dup_list)))]
dup_list |> kruskal.test() -> dup_test
dup_list |> dunns_test() -> dup.dunn.test
write.csv(as.data.frame(dup.dunn.test$test.result), "../../dicty_data/expr_dunns_test.csv")

# Enrich vs. dS human:
vec_list_human <- jsonlite::fromJSON(means_vecs_human_file)
names(vec_list_human$K9me3) <- as.character(1:10)
vec_list_human$K9me3 |> dunns_test() ->  human.dunn.test
medians <- lapply(vec_list_human$K9me3, median)
write.csv(as.data.frame(human.dunn.test$test.result), "../../dicty_data/h3k9me3_dunns_test_humans.csv")
