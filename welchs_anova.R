library(jsonlite)

make_df <- function(quantile_labels, means_vec){
  df <- data.frame(
    Enrichment=unlist(means_vec),
    Quant = factor(rep(quantile_labels, times=sapply(means_vec, length)))
  )
  return(df)
}

means_vecs_expr_file <- "./data/julia_serialized/means_vecs_expr.json"
means_vecs_ds_file <- "./data/julia_serialized/means_vecs_ds.json"
means_vecs_ds_human_file <- "./data/julia_serialized/means_vecs_ds_human.json"

means_vecs_expr <- fromJSON(means_vecs_expr_file)
means_vecs_ds <- fromJSON(means_vecs_ds_file)
means_vecs_ds_human <- fromJSON(means_vecs_ds_human_file)

quantile_labels <- paste("q", 1:10, sep="")

h3k4me3_means <- means_vecs_expr$K4me3
h3k27ac_means <- means_vecs_expr$K27ac
h3k9me3_means <- means_vecs_expr$K9me3
atac_means <- means_vecs_expr$ATAC

h3k4me3_df <- make_df(quantile_labels, h3k4me3_means)
h3k27ac_df <- make_df(quantile_labels, h3k27ac_means)
h3k9me3_df <- make_df(quantile_labels, h3k9me3_means)
atac_df <- make_df(quantile_labels, atac_means)

h3k4me3_results <- oneway.test(Enrichment ~ Quant, data=h3k4me3_df, var.equal=FALSE)
h3k27ac_results <- oneway.test(Enrichment ~ Quant, data=h3k27ac_df, var.equal=FALSE)
h3k9me3_results <- oneway.test(Enrichment ~ Quant, data=h3k9me3_df, var.equal=FALSE)
atac_results <- oneway.test(Enrichment ~ Quant, data=atac_df, var.equal=FALSE)

h3k9me3_means_human <- means_vecs_ds_human$K9me3
h3k9me3_df_human <- make_df(quantile_labels, h3k9me3_means_human)
h3k9me3_results_human <- oneway.test(Enrichment ~ Quant, data=h3k9me3_df_human, var.equal=FALSE)
