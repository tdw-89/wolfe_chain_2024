# Figure 4: low-dS paralogs (two lowest dS deciles) by overlap with H3K9me3, H3K27ac, H3K4me3 and ATAC peaks,
# with the expression of each intersection. Intersections that include H3K9me3 are highlighted; the
# horizontal line in each expression box marks the median of all genes with (or without) H3K9me3.
# Reads low_ds_df.csv (low_ds_paralogs.jl); writes ../../dicty_data/figures/fig4.pdf. The expression test
# is in low_ds_paralogs.jl (stats/low_ds_paralogs.csv).
suppressPackageStartupMessages({
  library(ComplexUpset)
  library(ggplot2)
})

data_dir <- "../../dicty_data"
fig_dir <- file.path(data_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)

low_ds <- read.csv(file.path(data_dir, "low_ds_df.csv"))
set_names <- c("H3K9me3", "H3K27ac", "H3K4me3", "ATAC")
names(low_ds) <- c("GeneID", set_names, "Expression")
low_ds[set_names] <- lapply(low_ds[set_names], function(x) as.logical(x))
low_ds$GroupMedian <- ave(low_ds$Expression, low_ds$H3K9me3, FUN = median)

k9_color <- "#9F313F"
plot <- upset(
  low_ds, set_names,
  name = "",
  base_annotations = list("Intersection size" = intersection_size(text = list(size = 3.5))),
  annotations = list(
    "Expression" = ggplot(mapping = aes(x = intersection, y = Expression)) +
      geom_boxplot(outlier.size = 0.8) +
      geom_errorbar(aes(ymin = GroupMedian, ymax = GroupMedian), color = "navyblue", linewidth = 0.8, width = 0.9) +
      ylab(expression(log[2] * "(TPM + 0.5)"))
  ),
  queries = list(upset_query(set = "H3K9me3", color = k9_color, fill = k9_color,
                             only_components = c("intersections_matrix", "Intersection size"))),
  themes = upset_modify_themes(list(default = theme(text = element_text(size = 12))))
)
ggsave(file.path(fig_dir, "fig4.pdf"), plot, width = 10, height = 7)
