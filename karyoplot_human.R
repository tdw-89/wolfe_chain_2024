library(karyoploteR)
library(rtracklayer)
library(ChIPseeker)
library(dplyr)
library(plyranges)
library(tools)
library(ggplot2)

# Files
# ensembl_99_gff <- "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.99.gff3"
ensembl_99_lengths <- "../../dicty_data/mammals/primates/h_sapiens/Ensembl_99/chromosome_lengths.txt"
repeat_data_file <- "../../dicty_data/mammals/primates/h_sapiens/GRCh38_UCSC_rmsk.tsv.gz"
peak_file <- "../../dicty_data/human_h3k9me3_avg.bed"

# Load the averaged peak data
peak_gr <- import.bed(peak_file)
# peak_gr <- peak_gr %>% filter(score >= 500)
peak_gr$score <- peak_gr$score / 500
seqlevels(peak_gr) <- paste("chr", seqlevels(peak_gr), sep="")

# Load the repeat data
repeat_data <- read.delim(gzfile(repeat_data_file), sep = "\t", header = T)
repeat_data <- repeat_data |> select(c(6, 7, 8, 12, 13))
names(repeat_data) <- c("Chromosome", "Start", "End", "Type", "Family")
#repeat_data <- repeat_data |> filter(!grepl("\\?", Type))
# repeat_data$Chromosome <- gsub("chr", "", repeat_data$Chromosome)
repeat_data <- repeat_data |> filter(Chromosome %in% seqlevels(peak_gr))

# Create sub-groups
rna_groups <- c("snRNA", "tRNA", "srpRNA", "rRNA", "scRNA", "RNA")
misc_groups <- c("Simple_repeat", "Satellite", "Low_complexity")
LINEs <- repeat_data |> filter(Type == "LINE")
L1s <- repeat_data |> filter(Family == "L1")
SINEs <- repeat_data |> filter(Type == "SINE")
Alus <- repeat_data |> filter(Family == "Alu")
active_TEs <- rbind(L1s, Alus)
LTRs <- repeat_data |> filter(Type == "LTR")
ERVs <- repeat_data |> filter(grepl("ERV", Family))
Retroposons <- repeat_data |> filter(Type == "Retroposon")
DNAs <- repeat_data |> filter(Type == "DNA")
retro_TEs <- rbind(LINEs, SINEs, LTRs, Retroposons)
all_TEs <- rbind(retro_TEs, DNAs)
RNAs <- repeat_data |> filter(Type %in% rna_groups)
other_repeats <-  repeat_data |> filter(Type %in% misc_groups)
unknowns <- repeat_data |> filter(Type == "Unknown")
dubious <- repeat_data |> filter(Type %in% c("DNA?", "LTR?", "SINE?", "RC?"))
rcs <- repeat_data |> filter(Type == "RC")

# Create the karyoplot
for(chr_nm in seqlevels(peak_gr)){
  min_start <- min(start(peak_gr[peak_gr@seqnames == chr_nm]))
  max_end <- max(end(peak_gr[peak_gr@seqnames == chr_nm]))
  total_span <- max_end - min_start + 1
  one_third <- ceiling(total_span / 3)
  zoom_start <- one_third
  zoom_end <- one_third + 500000
  kp <- plotKaryotype(genome = "hg38", 
                      plot.type = 2, 
                      # chromosomes = chr_nm),
                      zoom = toGRanges(chr_nm, zoom_start, zoom_end))
  kpPlotRegions(kp, peak_gr,
               col = "red", border = lighter("red"))
  # kpPlotDensity(kp, data = GRanges(L1s[(L1s$End - L1s$Start) >= 1500,]), col = "blue", border = "blue",
  #               data.panel = 2)
  kpPlotRegions(kp, data = GRanges(active_TEs[(active_TEs$End - active_TEs$Start) >= 1500,]), col = "blue", border = "blue",
                data.panel = 2)
}

