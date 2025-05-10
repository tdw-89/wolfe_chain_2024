library(dunn.test)

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
expr_deciles <- "/Users/thomasd.wolfe/Library/CloudStorage/GoogleDrive-eccwolfe@gmail.com/My Drive/UML/Grad_school/research/code/analysis_scripts/dicty/final_code/data/julia_serialized/expression_deciles.json"

# Enrich vs. dS:
vec_list <- jsonlite::fromJSON(means_vecs_file)
names(vec_list$K27ac) <- as.character(1:10)
names(vec_list$K4me3) <- as.character(1:10)
names(vec_list$K9me3) <- as.character(1:10)
names(vec_list$ATAC) <- as.character(1:10)

vec_list$K27ac |> kruskal.test() -> k27.kw.test
vec_list$K4me3 |> kruskal.test() -> k4.kw.test
vec_list$K9me3 |> kruskal.test() -> k9.kw.test
vec_list$ATAC |> kruskal.test() -> atac.kw.test
pvals_adj <- p.adjust(c(k27.kw.test$p.value, 
                        k4.kw.test$p.value, 
                        k9.kw.test$p.value, 
                        atac.kw.test$p.value), 
                      method = "bonferroni")

vec_list$K4me3 |> dunns_test() -> k4.dunn.test
write.csv(as.data.frame(k4.dunn.test$test.result), "./data/k4_dunns_test.csv")

# Expression vs. dS:
vec_list_expr <- jsonlite::fromJSON(expr_deciles)
dup_list <- vec_list_expr[-1]
dup_list <- dup_list[order(as.integer(names(dup_list)))]
dup_list |> kruskal.test() -> dup_test
dup_list |> dunns_test() -> dup.dunn.test

# Expression - Singletons vs dS deciles:
p1 <- wilcox.test(vec_list_expr$singleton, vec_list_expr$`1`)$p.value
p2 <- wilcox.test(vec_list_expr$singleton, vec_list_expr$`2`)$p.value
p3 <- wilcox.test(vec_list_expr$singleton, vec_list_expr$`3`)$p.value
p4 <- wilcox.test(vec_list_expr$singleton, vec_list_expr$`4`)$p.value
p5 <- wilcox.test(vec_list_expr$singleton, vec_list_expr$`5`)$p.value
p6 <- wilcox.test(vec_list_expr$singleton, vec_list_expr$`6`)$p.value
p7 <- wilcox.test(vec_list_expr$singleton, vec_list_expr$`7`)$p.value
p8 <- wilcox.test(vec_list_expr$singleton, vec_list_expr$`8`)$p.value
p9 <- wilcox.test(vec_list_expr$singleton, vec_list_expr$`9`)$p.value
p10 <- wilcox.test(vec_list_expr$singleton, vec_list_expr$`10`)$p.value

pvals <- c(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)
pvals_adj <- p.adjust(pvals, method = "BH")

fisher.test(matrix(c(6113, 756, 3685, 544), nrow=2, ncol=2))
fisher.test(matrix(c(4866, 149, 6113, 756), nrow=2, ncol=2))
fisher.test(matrix(c(2753, 170, 3469, 289), nrow=2, ncol=2))


