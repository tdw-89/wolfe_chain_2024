#!/usr/bin/env bash
# Run the analysis from the TE annotations to the statistics and figures, in dependency order.
#
# Statistics are written to ../../dicty_data/stats and manuscript figures to ../../dicty_data/figures.
# dS is from Ensembl Compara by default; analyses that are also reported with YN00 dS are run a second
# time with --yn00 (pairs and YN00 dS from old_pairs_with_new_ds_dicty.csv).
#
# Upstream inputs that are not rebuilt here: the paralog and singleton lists (paralog_filt.jl,
# human_paralog_filt.jl), filtered expression (expr_data_filt.jl), YN00 dS (ds_dist_comparison.qmd)
# the peak files, and the copy-number-adjusted expression (cnv_adj.jl, whose serialized input-DNA
# coverage was written by an older version of BioinfoTools and no longer deserializes).
set -euo pipefail
cd "$(dirname "$0")"

jl() { echo "== $*"; julia --project=. --threads=auto "$@"; }

# TE distance, TE density and the H3K9me3 / TE overlap test
jl te_distance.jl
jl te_distance_human.jl
jl h3k9me3_dist_dicty.jl
./te_windows.sh
./te_genome_fraction.sh
jl gene_counts.jl
jl te_density.jl
python3 h3k9me3_te_dist.py

# Chromatin features vs expression (writes filtered/final_gene_list.txt, read by later scripts)
for stage in All V M F; do jl enrichment_vs_expression.jl "$stage"; done

# Chromatin features vs dS
for stage in All V M F; do jl enrichment_vs_ds.jl "$stage"; done
jl enrichment_vs_ds.jl All --yn00

# Human: TE tests (also caches the H3K9me3 peak data) and H3K9me3 vs dS
jl enrichment_vs_te_dist_human.jl
jl enrichment_vs_ds_human.jl

# Expression vs dS, then Kruskal–Wallis / Dunn's tests on the decile vectors written above
jl expression_vs_ds.jl
jl expression_vs_ds.jl --yn00
Rscript dunns_test.R

# Low-dS paralogs, TE density tests, TE gene expression
jl low_ds_paralogs.jl
jl enrichment_vs_te_dist.jl
jl gene_type_expression.jl

# Regression of dS on chromatin features
for stage in All V M F; do jl mlr_dicty.jl "$stage"; done
jl mlr_dicty.jl All --yn00
for stage in All M F; do jl pca_dicty.jl "$stage"; done

# LaTeX rows for the manuscript tables (../../dicty_data/stats/latex)
jl latex_tables.jl

# Figures (PDF, in ../../dicty_data/figures) and the Ensembl vs YN00 report
jl make_figures.jl
Rscript low_ds_upset.R
Rscript karyoplot_dicty.R
quarto render ds_analysis_comparison.qmd
