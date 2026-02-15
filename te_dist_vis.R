library(karyoploteR)
library(circlize)
library(regioneR)
library(rtracklayer)
library(dplyr)
library(conflicted)
conflicts_prefer(dplyr::filter)

# Files
te_file <- "../../dicty_data/rm_genomes/ensembl_52/full/results_v1/onecodetofindthemall/Dictyostelium_discoideum.dicty_2.7.dna.toplevel_aggregated_compat.csv"
chrom_lengths_ensembl <- "../../dicty_data/AX4/genome_ver_2_7/ensembl_52/chromosome_lengths_ensembl.txt"
peak_dir <- "../../dicty_data/wang_et_al/processed/run_1_ensembl52"
peak_dir_2 <- "../../dicty_data/wang_et_al/processed/run_2_ensembl52"

# load dicty genome data:
chrom_length_df <- read.delim(chrom_lengths_ensembl, header=F)
chrom_length_df$V1 <- paste("chr", chrom_length_df$V1, sep="")
chrom_seqs <- Seqinfo(chrom_length_df$V1, seqlengths = chrom_length_df$V2)

# Load dicty transposon data:
te_df <- read.csv(te_file, header=T, stringsAsFactors = F)
te_df$Chromosome <- paste("chr", te_df$Chromosome, sep="")
dna_te_df <- te_df |> filter(startsWith(Type, "DNA"))
rte_df <- te_df |> filter(!startsWith(Type, "DNA"))
names(rte_df) <- c("chr", "start", "end", "type", "family", "id")
names(dna_te_df) <- c("chr", "start", "end", "type", "family", "id")

##### Create initial plot #################################################################################
chroms_to_show <- paste("chr", as.character(seq(1:6)), sep="")
kp <- plotKaryotype(genome = chrom_seqs, plot.type = 2, chromosomes = chroms_to_show) #####################
##### Create initial plot #################################################################################

window_size <- 10000

rte_grange <- toGRanges(rte_df)
dna_te_grange <- toGRanges(dna_te_df)

# Add to the plots
kpPlotDensity(kp, data = dna_te_grange, col="#AACCFF", window.size = window_size)
# kpPlotRegions(kp, data = dna_te_grange, avoid.overlapping=F, col="#AACCFF")

kpPlotDensity(kp, data = rte_grange, col="#EED0A0", window.size = window_size)
# kpPlotRegions(kp, data = rte_grange, avoid.overlapping=F, col="#EED0A0")

# Load K9me3 data
peak_files <- list.files(peak_dir, full.names = T)
k9me3_files <- peak_files[grepl("9me3", peak_files)]
k9me3_dfs <- lapply(k9me3_files, read.delim, header=F)
k9me3_dfs_trimmed <- lapply(k9me3_dfs, function(df){
  new_df <- df[,c(1,2,3)]
  names(new_df) <- c("chr","start","end")
  return(new_df)
})

names(k9me3_dfs_trimmed) <- basename(k9me3_files)

# Rename the chromosomes to the simple name (e.g., chrX)
for(d in 1:length(k9me3_dfs_trimmed)){
  
  for(i in 1:nrow(k9me3_dfs_trimmed[[d]])){
    temp_chr <- k9me3_dfs_trimmed[[d]]$chr[i]
    if(!is.na(as.numeric(temp_chr))){
      k9me3_dfs_trimmed[[d]]$chr[i] <- paste("chr", temp_chr, sep="")
    }
  }
}

# Convert dfs to individual Granges
k9me3_granges <- lapply(k9me3_dfs_trimmed, toGRanges)

# Convert to combined Grange
k9me3_df_combined <- bind_rows(k9me3_dfs_trimmed)
k9me3_grange_combined <- toGRanges(k9me3_df_combined)
# kpPlotRegions(kp, k9me3_grange_combined)

# Add to plot:
kpPlotDensity(kp, k9me3_grange_combined, col="#9F313F", window.size = window_size, data.panel = 2, r0=0, r1=0.5)
# kpPlotRegions(kp, k9me3_grange_combined, col="#9F313F", avoid.overlapping=F, data.panel = 2, r0=0, r1=0.5)

# Load ATAC data:
peak_files_2 <- list.files(peak_dir_2, full.names = T)
atac_files <- peak_files_2[grepl("1E6", peak_files_2)]
# Remove 'S' sample
atac_files <- atac_files[!(grepl("-SA-", atac_files))]
atac_dfs <- lapply(atac_files, read.delim, header=F)
atac_dfs_trimmed <- lapply(atac_dfs, function(df){
  new_df <- df[,c(1,2,3)]
  names(new_df) <- c("chr","start","end")
  return(new_df)
})

names(atac_dfs_trimmed) <- basename(atac_files)

# Rename the chromosomes to the simple name (e.g., chrX)
for(d in 1:length(atac_dfs_trimmed)){
  
  for(i in 1:nrow(atac_dfs_trimmed[[d]])){
    
    conv_ind <- grepl(atac_dfs_trimmed[[d]]$chr[i], chrom_conversion$ID)
    atac_dfs_trimmed[[d]]$chr[i] <- chrom_conversion$name[conv_ind]
  }
}

# Convert to Grange
atac_granges <- lapply(atac_dfs_trimmed, toGRanges)
names(atac_granges) <- names(atac_dfs_trimmed)
atac_df_combined <- bind_rows(atac_dfs_trimmed)
atac_grange_combined <- toGRanges(atac_df_combined)

# Add to plot:
kpPlotDensity(kp, atac_grange_combined, col="#006F1F", window.size = 30000, data.panel = 2, r0=0.6, r1=1)
# kpPlotRegions(kp, atac_grange_combined, col="#006F1F", avoid.overlapping=F, data.panel = 2, r0=0.6, r1=1)

# Gray-out the ax4 dup region in the kp plot
# NOTE: duplication info: http://dictybase.org/Downloads/
dup_region_df <- data.frame("chr" = c("chr2"), "start" = 2263132, "end" = 3768654)
dup_region_grange <- toGRanges(dup_region_df)
kpPlotRegions(kp, dup_region_grange, col="#808080", data.panel = 2)
kpPlotRegions(kp, dup_region_grange, col="#808080", data.panel = 1)

#### INDIVIDUAL CHROMS ####

for(chrom_grange in chrom_granges){
  
  if(chrom_grange@ranges@width[1] >= 30000){
    
    kp_ind <- plotKaryotype(genome = chrom_grange, plot.type = 2, main = "All tissue types combined")
    kpPlotRegions(kp_ind, data = rte_grange, avoid.overlapping=F, col="#EED0A0")
    kpPlotRegions(kp_ind, data = dna_te_grange, avoid.overlapping=F, col="#AACCFF")
    kpPlotRegions(kp_ind, k9me3_grange_combined, avoid.overlapping=F, col="#9F313F", data.panel = 2)
    kpPlotDensity(kp_ind, atac_grange_combined, window.size = 30000, col="#006F1F", data.panel = 2)
    kpPlotRegions(kp_ind, dup_region_grange, avoid.overlapping=F, col="#808080", data.panel = 2)
    kpPlotRegions(kp_ind, dup_region_grange, avoid.overlapping=F, col="#808080", data.panel = 1)
  }
}

# Individual cell types:
for(chrom_grange in chrom_granges){
  
  if(chrom_grange@ranges@width[1] >= 30000){
    
    kp_ind <- plotKaryotype(genome = chrom_grange, plot.type = 2, main = "F")
    kpPlotRegions(kp_ind, data = rte_grange, avoid.overlapping=F, col="#EED0A0")
    kpPlotRegions(kp_ind, data = dna_te_grange, avoid.overlapping=F, col="#AACCFF")
    kpPlotRegions(kp_ind, k9me3_granges$FA_K9me3_fastq_dup_removed_converted.narrowPeak, avoid.overlapping=F, col="#9F313F", data.panel = 2)
    kpPlotDensity(kp_ind, atac_granges$`171031-13-FA-1E6.narrowPeak`, window.size = 30000, col="#006F1F", data.panel = 2)
    kpPlotRegions(kp_ind, dup_region_grange, avoid.overlapping=F, col="#808080", data.panel = 2)
    kpPlotRegions(kp_ind, dup_region_grange, avoid.overlapping=F, col="#808080", data.panel = 1)
  }
}

for(chrom_grange in chrom_granges){
  
  if(chrom_grange@ranges@width[1] >= 30000){
    
    kp_ind <- plotKaryotype(genome = chrom_grange, plot.type = 2, main = "M")
    kpPlotRegions(kp_ind, data = rte_grange, avoid.overlapping=F, col="#EED0A0")
    kpPlotRegions(kp_ind, data = dna_te_grange, avoid.overlapping=F, col="#AACCFF")
    kpPlotRegions(kp_ind, k9me3_granges$MA_K9me3_fastq_dup_removed_converted.narrowPeak, avoid.overlapping=F, col="#9F313F", data.panel = 2)
    kpPlotDensity(kp_ind, atac_granges$`171031-11-MA-1E6.narrowPeak`, window.size = 30000, col="#006F1F", data.panel = 2)
    kpPlotRegions(kp_ind, dup_region_grange, avoid.overlapping=F, col="#808080", data.panel = 2)
    kpPlotRegions(kp_ind, dup_region_grange, avoid.overlapping=F, col="#808080", data.panel = 1)
  }
}

for(chrom_grange in chrom_granges){
  
  if(chrom_grange@ranges@width[1] >= 30000){
    
    kp_ind <- plotKaryotype(genome = chrom_grange, plot.type = 2, main = "V")
    kpPlotRegions(kp_ind, data = rte_grange, avoid.overlapping=F, col="#EED0A0")
    kpPlotRegions(kp_ind, data = dna_te_grange, avoid.overlapping=F, col="#AACCFF")
    kpPlotRegions(kp_ind, k9me3_granges$VA_K9me3_fastq_dup_removed_converted.narrowPeak, avoid.overlapping=F, col="#9F313F", data.panel = 2)
    kpPlotDensity(kp_ind, atac_granges$`171031-7-VA-1E6_S10.narrowPeak`, window.size = 30000, col="#006F1F", data.panel = 2)
    kpPlotRegions(kp_ind, dup_region_grange, avoid.overlapping=F, col="#808080", data.panel = 2)
    kpPlotRegions(kp_ind, dup_region_grange, avoid.overlapping=F, col="#808080", data.panel = 1)
  }
}

