library(rtracklayer)
library(plyranges)


duplicate_file <- "./data/filtered/human_paralog_info_filt.csv"
blacklist_file <- "../../../../data/mammals/primates/h_sapiens/ENCODE_histone_mods/ENCFF356LFX.bed.gz"
grch38_ensembl_file <- "../../../../data/mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.99.gff3"
repeat_data_file <- "../../../../data/mammals/primates/h_sapiens/GRCh38_UCSC_rmsk.tsv.gz"


perc_overlap <- function(dup_start, dup_end, black_start, black_end){
  overlap_len <- length(intersect(dup_start:dup_end, black_start:black_end))
  gene_len <- dup_end - dup_start + 1
  return(overlap_len/gene_len * 100)
}

len_overlap <- function(dup_start, dup_end, black_start, black_end){
  overlap_len <- length(intersect(dup_start:dup_end, black_start:black_end))
  return(overlap_len)
}

######################
# DUPLICATE OVERLAPS #
######################

dup_df <- read.csv(duplicate_file)
blacklist_gr <- import.bed(blacklist_file)
grch38_gr <- import.gff(grch38_ensembl_file)

dup_ids <- c(dup_df$GeneID, dup_df$ParalogID)
dup_gr <- grch38_gr |> filter(gene_id %in% dup_ids)
seqlevels(dup_gr) <- paste("chr", seqlevels(dup_gr), sep="")

overlaps <- findOverlaps(dup_gr, blacklist_gr)
from_df <- as.data.frame(dup_gr[queryHits(overlaps)])
names(from_df) <- paste(names(from_df), "_dup", sep="")
to_df <- as.data.frame(blacklist_gr[subjectHits(overlaps)])

overlap_df <- cbind(from_df, to_df)
overlap_df$perc_dup_overlapped <- mapply(
  perc_overlap,
  overlap_df$start_dup,
  overlap_df$end_dup,
  overlap_df$start,
  overlap_df$end
)

###################
# REPEAT OVERLAPS #
###################

repeat_data_df <- read.delim(gzfile(repeat_data_file), sep = "\t", header = T)
repeat_data_df <- repeat_data_df |> select(c(6, 7, 8, 12, 13))
names(repeat_data_df) <- c("Chromosome", "Start", "End", "Type", "Family")
repeat_data_df <- repeat_data_df |> filter(Chromosome %in% seqlevels(blacklist_gr))

repeat_types <- unique(repeat_data_df$Type)
overlap_df <- data.frame("RepeatType" = repeat_types, 
                         "NumOverlaps" = -1, 
                         "NumRepeats" = -1, 
                         "PercRepeatsOverlapped" = -1,
                         "NumNucleotidesOverlapped" = -1,
                         "NumNucleotides" = -1,
                         "PercNucleotidesOverlapped" = -1)
for(i in 1:length(repeat_types)){
  print(i)
  repeat_type <- repeat_types[i]
  repeat_gr <- repeat_data_df |> filter(Type == repeat_type)
  repeat_gr <- makeGRangesFromDataFrame(repeat_gr, keep.extra.columns = T)
  overlaps <- findOverlaps(repeat_gr, blacklist_gr)
  from_df <- as.data.frame(repeat_gr[queryHits(overlaps)])
  names(from_df) <- paste(names(from_df), "_repeat", sep="")
  to_df <- as.data.frame(blacklist_gr[subjectHits(overlaps)])
  
  if(isEmpty(from_df) || isEmpty(to_df)){
    overlap_df$NumOverlaps[i] <- 0
    overlap_df$NumRepeats[i] <- length(repeat_gr)
    overlap_df$PercRepeatsOverlapped[i] <- 0
    overlap_df$NumNucleotides[i] <- sum(width(repeat_gr))
    overlap_df$NumNucleotidesOverlapped[i] <- 0
    overlap_df$PercNucleotidesOverlapped[i] <- 0
    print(paste(i, "done"))
    next
  }else{
    overlap_df$NumOverlaps[i] <- length(unique(queryHits(overlaps)))
    overlap_df$NumRepeats[i] <- length(repeat_gr)
    overlap_df$PercRepeatsOverlapped[i] <- overlap_df$NumOverlaps[i] / overlap_df$NumRepeats[i] * 100
    overlap_df$NumNucleotides[i] <- sum(width(repeat_gr))
    overlap_df$NumNucleotidesOverlapped[i] <- sum(mapply(
      len_overlap,
      from_df$start_repeat,
      from_df$end_repeat,
      to_df$start,
      to_df$end
    ))
    overlap_df$PercNucleotidesOverlapped[i] <- overlap_df$NumNucleotidesOverlapped[i] / overlap_df$NumNucleotides[i] * 100
  }
  
  print(paste(i, "done"))
}

write.csv(overlap_df, "./data/te_blacklist_overlap_info.csv", row.names = F)


