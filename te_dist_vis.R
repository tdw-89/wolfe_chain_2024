library(karyoploteR)
library(circlize)
library(regioneR)
library(rtracklayer)
library(dplyr)
conflicts_prefer(dplyr::filter)

# Files
te_files_dir <- "../../dicty_data/AX4/genome_ver_2_7/TEs"
gff_dir <- "../../dicty_data/AX4/genome_ver_2_7/gff3"
peak_dir <- "../../dicty_data/wang_et_al/processed/run_1_ensembl52"
peak_dir_2 <- "../../dicty_data/wang_et_al/processed/run_2_ensembl52"

# Functions: #
trim_names <- function(x){
  temp_name <- strsplit(x, "_")[[1]][2]
  return(paste("chr", temp_name, sep=""))
}

# load dicty genome data:
dicty_gff_files <- list.files(gff_dir, pattern = "chromosome.*gff", full.names = T)
dicty_granges <- lapply(dicty_gff_files, rtracklayer::readGFFAsGRanges)
names(dicty_granges) <- unlist(lapply(dicty_gff_files, function(x){return(strsplit(basename(x), "[.]")[[1]][1])}))
dicty_granges_list <- GRangesList(dicty_granges)

chrom_df <- data.frame("chr" = unname(sapply(names(dicty_granges), trim_names)),
                       "start" = rep(1, length(dicty_granges)),
                       "end" = unname(sapply(dicty_granges, function(x){x@ranges@width[1]})))

chrom_ids <- levels(dicty_granges_list@unlistData@seqnames)
chrom_conversion <- data.frame("ID" = chrom_ids, "name" = chrom_df$chr)

chrom_df <- chrom_df %>% filter(chr %in% c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6"))
genome_grange <- toGRanges(chrom_df)
chrom_granges <- lapply(chrom_df$chr, function(x, df){
  temp_df <- df[df$chr == x,]
  return(toGRanges(temp_df))
}, chrom_df)

names(chrom_granges) <- chrom_df$chr

# Load dicty transposon data:
te_bed_files <- list.files(te_files_dir, pattern = "ranges.bed", full.names = T)
rte_df <- read.delim(te_bed_files[grepl("rte", te_bed_files)], header=F)
dna_te_df <- read.delim(te_bed_files[grepl("dna_", te_bed_files)], header=F)
names(rte_df) <- c("chr", "start", "end", "name", "score", "strand")
names(dna_te_df) <- c("chr", "start", "end", "name", "score", "strand")

# Rename the chromosomes to the simple name (e.g., chrX)

for(i in 1:nrow(rte_df)){
  
  conv_ind <- grepl(rte_df$chr[i], chrom_conversion$ID)
  rte_df$chr[i] <- chrom_conversion$name[conv_ind]
}
for(i in 1:nrow(dna_te_df)){
  
  conv_ind <- grepl(dna_te_df$chr[i], chrom_conversion$ID)
  dna_te_df$chr[i] <- chrom_conversion$name[conv_ind]
}


##### Create initial plot ######################################################
kp <- plotKaryotype(genome = genome_grange, plot.type = 2) #####################
##### Create initial plot ######################################################

rte_grange <- toGRanges(rte_df)
dna_te_grange <- toGRanges(dna_te_df)

# Add to the plots
kpPlotRegions(kp, data = rte_grange, avoid.overlapping=F, col="#EED0A0")
# kpPlotDensity(kp, data = rte_grange)

kpPlotRegions(kp, data = dna_te_grange, avoid.overlapping=F, col="#AACCFF")
# kpPlotDensity(kp_dna_te, data = dna_te_grange)

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

k9me3_dfs_trimmed <- lapply(k9me3_dfs_trimmed, function(df){
  df |> filter(chr %in% chrom_conversion$name)
})

# Convert dfs to individual Granges
k9me3_granges <- lapply(k9me3_dfs_trimmed, toGRanges)

# Convert to combined Grange
k9me3_df_combined <- bind_rows(k9me3_dfs_trimmed)
k9me3_grange_combined <- toGRanges(k9me3_df_combined)
# kpPlotRegions(kp, k9me3_grange_combined)

# Add to plot:
kpPlotRegions(kp, k9me3_grange_combined, col="#9F313F", avoid.overlapping=F, data.panel = 2, r0=0, r1=0.5)

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

