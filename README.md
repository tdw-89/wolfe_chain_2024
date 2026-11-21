# Epigenetic components of paralog evolution in *D. discoideum*

Analysis code for the *D. discoideum* / human paralog manuscript. Data are read from and written to
`../../dicty_data` (relative to this directory).

## Running

```sh
./run_analysis.sh
```

runs everything from the TE annotations to the statistics and figures, in dependency order. It needs
Julia (the project in this directory), R with `dunn.test`, `jsonlite`, `ComplexUpset`, `ggplot2` and
`karyoploteR`, Python 3, bedtools (`$BEDTOOLS`, then `PATH`, then `~/Documents/bedtools2/bin`) and Quarto.

Scripts can also be run one at a time with `julia --project=. <script>.jl [stage] [--yn00]`:

- `stage`: life-cycle stage, `All` (default), `V`, `M` or `F`. `All` averages the vegetative, mound and
  fruiting-body samples (ChIP-seq) and the streaming, mound and fruiting-body samples (ATAC-seq).
- `--yn00`: use YN00 dS (`old_pairs_with_new_ds_dicty.csv`) instead of Ensembl Compara dS, for the
  comparison of the two estimates. Output files get a `_yn00` (stats) or `_newds` / `_new_ds` (intermediate
  files) suffix. The manuscript uses Ensembl Compara dS.

Permutation tests use a fixed seed (set in `prelude.jl`).

## Statistics

Every statistic reported in the manuscript is written to `../../dicty_data/stats/`:

| File | Contents | Produced by |
|---|---|---|
| `gene_counts.csv`, `genes_overlapping_te.csv` | Gene counts and genes overlapping a TE (Table 1) | `gene_counts.jl` |
| `te_genome_fraction.csv` | Non-redundant TE bases by class (Table 1) | `te_genome_fraction.sh` |
| `te_coverage_human_annotations.csv` | Summed human TE annotation lengths (overlaps counted twice) | `te_distance_human.jl` |
| `h3k9me3_te_overlap.csv` | H3K9me3 / TE overlap and permutation test | `h3k9me3_te_dist.py` |
| `te_gene_expression.csv` | TE vs non-TE gene expression (Fig S2) | `gene_type_expression.jl` |
| `enrichment_vs_expression_<stage>.csv` | Correlations and Kruskal–Wallis tests, expression (Table 2) | `enrichment_vs_expression.jl` |
| `enrichment_vs_ds_<stage>[_yn00].csv` | Correlations and Kruskal–Wallis tests, dS (Table 2) | `enrichment_vs_ds.jl` |
| `enrichment_vs_ds_human.csv` | Human H3K9me3 vs dS | `enrichment_vs_ds_human.jl` |
| `kruskal_wallis_by_ds.csv`, `*_dunn.csv` | Kruskal–Wallis and Dunn's tests by dS decile (Fig 3, Table S1) | `dunns_test.R` |
| `expression_vs_ds*[_yn00].csv` | Expression vs dS: tests (including the trend after decile 1), decile summaries, paralogs vs singletons (Fig 3) | `expression_vs_ds.jl` |
| `low_ds_paralogs*.csv` | Low-dS paralogs with vs without H3K9me3 (Fig 4) | `low_ds_paralogs.jl` |
| `te_density_tests.csv` | TE density by paralog and H3K9me3 status | `enrichment_vs_te_dist.jl` |
| `te_tests_human_TE.csv` | Human TE distance, overlap and triplosensitivity tests (Table 4) | `enrichment_vs_te_dist_human.jl` |
| `mlr_<stage>_<ensembl\|yn00>[_model\|_vif].csv` | Regression of dS (Table 3, VIF table, per-stage tables) | `mlr_dicty.jl` |
| `pca_regression_<stage>[_model].csv` | Principal component regression | `pca_dicty.jl` |

`latex_tables.jl` turns these into LaTeX table rows (Tables 2, 3, S1 and the VIF, per-stage regression and
PC regression tables) in `../../dicty_data/stats/latex/`.

## Figures

Manuscript figures are written as vector PDF to `../../dicty_data/figures/`, named as in the manuscript
(`fig1.pdf` … `figS17.pdf`):

| Figures | Produced by |
|---|---|
| 1, 2, 3, 5, S3–S5, S7–S9, S13 | `make_figures.jl` |
| 4 | `low_ds_upset.R` |
| S1 | `karyoplot_dicty.R` |
| S2 | `gene_type_expression.jl` |
| S6, S11 | `enrichment_vs_te_dist.jl` |
| S10 | `expression_vs_ds.jl` |
| S12 | `enrichment_vs_te_dist_human.jl` |
| S14 | `mlr_dicty.jl` |
| S15–S17 | `pca_dicty.jl` |

The shared style (fonts, colours, sizes) is in `figure_style.jl`. Exploratory HTML plots are also written
to `../../dicty_data/` and `../../dicty_data/saved_figures/`.

## Inputs not rebuilt by `run_analysis.sh`

Paralog and singleton lists (`paralog_filt.jl`, `human_paralog_filt.jl`), filtered expression
(`expr_data_filt.jl`), YN00 dS (`ds_dist_comparison.qmd`), peak calls, and the copy-number-adjusted
expression from `cnv_adj.jl` (its serialized input-DNA coverage was written by an older BioinfoTools and no
longer deserializes; `filtered/expr_data_filt_adj.tsv` is used as is).
