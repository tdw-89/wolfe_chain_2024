library(dplyr)
library(ComplexUpset)
library(ggplot2)

# Function to calculate averages
calc_averages <- function(df, sets, val_col) {
  comb <- expand.grid(rep(list(0:1), length(sets)))
  colnames(comb) <- sets
  comb <- comb[rowSums(comb) > 0, ]
  comb$avg_value <- apply(comb, 1, function(row) {
    subset <- df[apply(df[sets] == row, 1, all), ]
    mean(subset[, val_col], na.rm = TRUE)
  })
  comb
}


# low_ds_table <- "./data/all_ds_df.csv"
low_ds_table <- "./data/low_ds_df.csv"

# Load the data as a df
low_ds_df <- read.csv(low_ds_table)
low_ds_df[, 2:5] <- lapply(low_ds_df[, 2:5], function(col) {
  return(as.integer(factor(col)) - 1)
})

set_names <- c("HasK9me3",
               "HasK27ac",
               "HasK4me3",
               "HasATAC")

expr_avgs <- calc_averages(low_ds_df, set_names, 6)

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
        "Expression" = ggplot(mapping=aes(x = intersection, y = Expr)) +
          geom_boxplot(lwd=1.08)
      ),
      themes = plt_theme)
