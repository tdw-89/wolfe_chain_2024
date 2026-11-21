# Figure S1: density of class I and class II TEs (above each chromosome) and of H3K9me3 peaks in any of the
# vegetative, mound and fruiting-body samples (below) along the six D. discoideum chromosomes.
# Writes ../../dicty_data/figures/figS01.pdf.
suppressPackageStartupMessages({
  library(karyoploteR)
  library(regioneR)
})

data_dir <- "../../dicty_data"
te_file <- file.path(data_dir, "rm_genomes/ensembl_52/full/results_v3/onecodetofindthemall",
                     "Dictyostelium_discoideum.dicty_2.7.dna.toplevel_aggregated_compat.csv")
chrom_lengths_file <- file.path(data_dir, "AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt")
peak_dir <- file.path(data_dir, "wang_et_al/processed/run_1_ensembl52")
fig_dir <- file.path(data_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)

window_size <- 10000
chroms <- paste0("chr", 1:6)

chrom_lengths <- read.delim(chrom_lengths_file, header = FALSE)
genome <- toGRanges(data.frame(chr = paste0("chr", chrom_lengths$V1), start = 1, end = chrom_lengths$V2))

to_granges <- function(chrom, start, end) toGRanges(data.frame(chr = paste0("chr", chrom), start = start, end = end))

tes <- read.csv(te_file, stringsAsFactors = FALSE)
is_dna <- startsWith(tes$Type, "DNA")
dna_tes <- to_granges(tes$Chromosome[is_dna], tes$Start[is_dna], tes$End[is_dna])
retro_tes <- to_granges(tes$Chromosome[!is_dna], tes$Start[!is_dna], tes$End[!is_dna])

k9_files <- list.files(peak_dir, pattern = "_[VMF]A_K9me3_.*\\.narrowPeak$", full.names = TRUE)
k9_peaks <- do.call(rbind, lapply(k9_files, read.delim, header = FALSE))
k9 <- to_granges(k9_peaks$V1, k9_peaks$V2, k9_peaks$V3)

pdf(file.path(fig_dir, "figS01.pdf"), width = 10, height = 6.5)
kp <- plotKaryotype(genome = genome, plot.type = 2, chromosomes = chroms)
kpPlotDensity(kp, data = dna_tes, col = "#AACCFF", window.size = window_size)
kpPlotDensity(kp, data = retro_tes, col = "#EED0A0", window.size = window_size)
kpPlotDensity(kp, data = k9, col = "#9F313F", window.size = window_size, data.panel = 2)
invisible(dev.off())
