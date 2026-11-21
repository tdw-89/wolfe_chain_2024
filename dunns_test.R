# Kruskal–Wallis and Dunn's post-hoc tests (Holm-adjusted, two-sided) on the decile vectors written by
# enrichment_vs_ds.jl, expression_vs_ds.jl and enrichment_vs_ds_human.jl. Results go to dicty_data/stats.
# Each analysis is run with Ensembl Compara dS and, where the YN00 files exist, with YN00 dS.
library(dunn.test)

ser_dir <- "../../dicty_data/julia_serialized"
stats_dir <- "../../dicty_data/stats"
dir.create(stats_dir, showWarnings = FALSE)

# groups: named list of numeric vectors (names are decile numbers)
run_tests <- function(groups, name) {
  kw <- kruskal.test(groups)
  # altp = TRUE gives two-sided p-values (the default reports one-sided p, to be compared with alpha/2)
  invisible(capture.output(d <- dunn.test(groups, method = "holm", altp = TRUE, kw = FALSE, table = FALSE, list = FALSE)))
  # dunn.test labels groups by position; map them back to the decile names
  labels <- sapply(strsplit(d$comparisons, " - "), function(pair) {
    paste(names(groups)[as.integer(pair)], collapse = " - ")
  })
  write.csv(data.frame(comparison = labels, Z = d$Z, p = d$altP, p_adj = d$altP.adjusted),
            file.path(stats_dir, paste0(name, "_dunn.csv")), row.names = FALSE)
  data.frame(test = name, chi_squared = unname(kw$statistic), df = unname(kw$parameter),
             p = kw$p.value, n = sum(lengths(groups)),
             medians = paste(signif(sapply(groups, median), 3), collapse = ";"))
}

deciles <- as.character(1:10)
kw_results <- list()
for (variant in c("", "_yn00")) {
  # H3K4me3 and ATAC coverage by dS decile, the marks with a significant Kruskal-Wallis test (Table S1)
  ds_file <- file.path(ser_dir, paste0("means_vecs_ds_All", if (variant == "") "" else "_newds", ".json"))
  if (file.exists(ds_file)) {
    ds_vecs <- jsonlite::fromJSON(ds_file)
    for (mark in c("K4me3", "ATAC")) {
      groups <- ds_vecs[[mark]]
      names(groups) <- deciles
      name <- paste0(c(K4me3 = "h3k4me3", ATAC = "atac")[[mark]], "_by_ds", variant)
      kw_results[[length(kw_results) + 1]] <- run_tests(groups, name)
    }
  }

  # Pair-average expression by dS decile (Figure 3); the file also holds the singletons, which are
  # not part of this test
  expr_file <- file.path(ser_dir, paste0("expression_deciles", if (variant == "") "" else "_new_ds", ".json"))
  if (file.exists(expr_file)) {
    expr <- jsonlite::fromJSON(expr_file)[deciles]
    kw_results[[length(kw_results) + 1]] <- run_tests(expr, paste0("expression_by_ds", variant))
  }

  # Human H3K9me3 coverage by dS decile (Table S1)
  human_file <- file.path(ser_dir, paste0("means_vecs_dS_human", if (variant == "") "" else "_newds", ".json"))
  if (file.exists(human_file)) {
    k9 <- jsonlite::fromJSON(human_file)$K9me3
    names(k9) <- deciles
    kw_results[[length(kw_results) + 1]] <- run_tests(k9, paste0("h3k9me3_by_ds_human", variant))
  }
}
write.csv(do.call(rbind, kw_results), file.path(stats_dir, "kruskal_wallis_by_ds.csv"), row.names = FALSE)
