library(rtracklayer)
library(tidyverse)

overlap_prop <- function(start1, end1, start2, end2){
  total_len <- max(end1, end2) - min(start1, start2) + 1
  intersect_len <- length(intersect(start1:end1, start2:end2))
  return(intersect_len/total_len)
}

# Files
sens_data_file <- "../../../../data/mammals/primates/h_sapiens/NIHMS1819123-supplement-14-dosage-sensitivity-predictions.csv"
lift_over_chain_file <- "../../../../data/mammals/primates/h_sapiens/GRCh37_to_GRCh38_Ensembl.chain"
filtered_paralog_data <- "../../../../code/analysis_scripts/dicty/final_code/data/filtered/human_paralog_info_filt.csv"
ensembl_99_gff <- "../../../../data/mammals/primates/h_sapiens/Ensembl_99/Homo_sapiens.GRCh38.99.gff3"
gencode_v19_gtf <- "../../../../data/mammals/primates/h_sapiens/gencode.v19.chr_patch_hapl_scaff.annotation.gtf"

# Load the GENCODE v19 GTF data
gtf_df_gene <- as.data.frame(rtracklayer::import(gencode_v19_gtf)) %>%
  filter(type == "gene" & gene_type == "protein_coding")

# Load the sensitivity data
sens_df <- read.csv(sens_data_file)

# Make sure all gene names are in the GTF data
isEmpty(setdiff(sens_df$Gene, unique(gtf_df_gene$gene_name)))

# Match gene symbols to coordinates
sens_df_j <- sens_df %>%
  full_join(gtf_df_gene, by = c(Gene = "gene_name"))

# Filter to rows with sensitivity data
sens_df_j <- sens_df_j %>%
  filter(!is.na(pTriplo))

# To handle duplicate entries, first group by gene name, and then take the
# entry (possibly multiple) with the lowest 'level' among the duplicates.
sens_df_j_no_dup <- sens_df_j %>%
  group_by(Gene) %>%
  filter(level == min(level)) %>%
  slice(1)

# The dosage sensitivity study only kept genes on autosomes:
num_chroms <- paste("chr", as.character(1:22), sep="")
sens_df_j_no_dup <- sens_df_j_no_dup[sens_df_j_no_dup$seqnames %in% num_chroms,]
sens_df_j_no_dup$seqnames <- sapply(sens_df_j_no_dup$seqnames, 
                                    function(str){return(gsub("chr", "", str))})

# Load the filtered paralog data, and identify the genomic coordinates of each
# paralog using the gff file
gff_df_99 <- as.data.frame(rtracklayer::import(ensembl_99_gff)) %>%
  filter(type == "gene" & biotype == "protein_coding")
dup_genes <- read.csv(filtered_paralog_data)
dup_df <- gff_df_99 %>%
  filter(gene_id %in% dup_genes$GeneID | gene_id %in% dup_genes$ParalogID)

# Create a genomic ranges object from the sensitivity/GTF combined dataset
sens_df_gr <- GRanges(
  seqnames = sens_df_j_no_dup$seqnames, ranges = IRanges(
    start = sens_df_j_no_dup$start,
    end = sens_df_j_no_dup$end
  ), strand = sens_df_j_no_dup$strand, gene_name = sens_df_j_no_dup$Gene,
  gene_id = sens_df_j_no_dup$gene_id, pTriplo = sens_df_j_no_dup$pTriplo
)

# Lift over the genomic ranges object to GRCh38
chain <- import.chain(lift_over_chain_file)
sens_df_gr_lifted <- unlist(liftOver(sens_df_gr, chain))


# Creat a genomic ranges object from the duplicate gene data
dup_gr <- GRanges(seqnames = dup_df$seqnames, ranges = IRanges(
  start = dup_df$start,
  end = dup_df$end
), strand = dup_df$strand, gene_name = dup_df$Name, gene_id = dup_df$gene_id)

# Find overlaps between the duplicate gene data and the lifted sensitivity data
overlaps <- findOverlaps(dup_gr, sens_df_gr_lifted)
df1 <- as.data.frame(dup_gr[queryHits(overlaps)])
names(df1) <- paste(names(df1), "_GRCh38", sep = "")
df2 <- as.data.frame(sens_df_gr_lifted[subjectHits(overlaps)])
names(df2) <- paste(names(df2), "_liftover_GRCh37", sep = "")
original_coords = sens_df_j_no_dup[match(df2$gene_id_liftover_GRCh37, 
                                         sens_df_j_no_dup$gene_id),] %>%
  ungroup %>%
  select(c(seqnames, start, end, width))
names(original_coords) <- paste(names(original_coords), "_original_GRCh37", sep = "")
overlap_df <- cbind(df1, df2, original_coords)
overlap_df$overlap_prop <- mapply(overlap_prop, 
                                  overlap_df$start_GRCh38, 
                                  overlap_df$end_GRCh38, 
                                  overlap_df$start_liftover_GRCh37, 
                                  overlap_df$end_liftover_GRCh37)

# Save the overlap df
write.csv(overlap_df, "./data/lift_over_data_full.csv")

# Filter to only the overlaps where >= 50% of the GRCh37 sequence overlaps the
# GRCh38 sequence
overlap_df$gene_id_liftover_GRCh37 <- unname(sapply(overlap_df$gene_id_liftover_GRCh37, function(str){return(strsplit(str, "\\.")[[1]][1])}))
overlap_df_filtered <- overlap_df[which(overlap_df$gene_id_GRCh38 == overlap_df$gene_id_liftover_GRCh37),]
overlap_df_filtered <- overlap_df_filtered %>% filter(overlap_prop > 0.5)
mismatches <- anti_join(overlap_df, overlap_df_filtered)
mismatches <- anti_join(mismatches, filter(mismatches, gene_id_GRCh38 %in% overlap_df_filtered$gene_id_GRCh38))

# Save the filtered overlaps and the mismatches
write.csv(overlap_df_filtered, "./data/lift_over_data_filtered.csv")
write.csv(mismatches, "./data/lift_over_data_mismatches.csv")

missing_ids <- dup_df %>% filter(!gene_id %in% overlap_df$gene_id_GRCh38)
