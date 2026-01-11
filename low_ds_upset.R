library(dplyr)
library(ComplexUpset)
library(ggplot2)
library(coro)

# low_ds_table <- "../../dicty_data/all_ds_df.csv"
low_ds_table <- "../../dicty_data/low_ds_df.csv"

# Load the data as a df
low_ds_df <- read.csv(low_ds_table)
set_names <- c("H3K9me3",
               "H3K27ac",
               "H3K4me3",
               "ATAC")
names(low_ds_df) <- c("Gene ID", set_names, "Expression")
low_ds_df[, 2:5] <- lapply(low_ds_df[, 2:5], function(col) {
  return(as.integer(factor(col)) - 1)
})

h3k9me3_med <- median(low_ds_df$Expression[low_ds_df$H3K9me3 == T])
non_h3k9me3_med <- median(low_ds_df$Expression[low_ds_df$H3K9me3 != T])
wilcox.test(low_ds_df$Expression[low_ds_df$H3K9me3 == T], 
            low_ds_df$Expression[low_ds_df$H3K9me3 != T])

medExpr <- generator(function() {
  for (x in c(
    non_h3k9me3_med,
    non_h3k9me3_med,
    non_h3k9me3_med,
    non_h3k9me3_med,
    non_h3k9me3_med,
    non_h3k9me3_med,
    h3k9me3_med,
    h3k9me3_med,
    h3k9me3_med,
    h3k9me3_med,
    h3k9me3_med,
    h3k9me3_med,
    non_h3k9me3_med
  )) {
    yield(x)
  }
})
f <- medExpr()

# Create the upset plot theme
    # Text
plt_theme <- upset_default_themes()
plt_theme$intersections_matrix[[1]]$text$size <- 18
plt_theme$`Intersection size`[[1]]$text$size <- 18
plt_theme$`Intersection size`[[1]]$axis.text$size <- 18
plt_theme$`Intersection size`[[1]]$axis.text.x$size <- 18
plt_theme$`Intersection size`[[1]]$axis.text.y$size <- 18
plt_theme$overall_sizes[[1]]$text$size <- 18
plt_theme$default[[1]]$text$size <- 18
  
upset(low_ds_df,
      set_names,
      annotations = list(
        "Expression" = ggplot(
            mapping=aes(
              x = intersection, 
              y = Expression
            )
          ) +
          geom_boxplot(lwd=1.08) +
          stat_summary(
            fun = NULL,
            fun.min = NULL,
            fun.max = f,
            geom = "errorbar", # Use errorbar to draw a horizontal line segment
            color = "navyblue",
            linewidth = 1.3,
            width = 1
          )
      ),
      themes = plt_theme)