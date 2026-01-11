library(conflicted)
library(MASS)
library(tidyverse)
library(pscl)

te_dist_file = "../../dicty_data/te_distance_human_TE.csv"
paralog_file = "../../dicty_data/filtered/human_paralog_info_filt.csv"
singleton_file = "../../dicty_data/filtered/human_singletons_filt.csv"
dosage_sens_file = "../../dicty_data/lift_over_data_filtered.csv"

# Load the distance data
dist_df <- read.csv(te_dist_file)
dist_df$TE_length <- dist_df$TE_end - dist_df$TE_start + 1
dist_df <- dplyr::select(dist_df, c(GeneID, Distance, TE_length))
summary(dist_df)

# Load the paralog and singleton data
paralog_df <- read.csv(paralog_file)
singleton_df <- read.csv(singleton_file)

# Load the dosage sensitivity data
dosage_sens_df <- read.csv(dosage_sens_file)

# singleton_dists = dist_df %>% dplyr::filter(GeneID %in% singleton_df$GeneID) %>% dplyr::filter(GeneID %in% dosage_sens_df$gene_id_GRCh38)
duplicate_dists <- dist_df %>% dplyr::filter(GeneID %in% c(paralog_df$GeneID, paralog_df$ParalogID))
duplicate_dosage <- dosage_sens_df %>% dplyr::filter(gene_id_GRCh38 %in% c(paralog_df$GeneID, paralog_df$ParalogID))
duplicate_dosage <- duplicate_dosage %>% dplyr::select(c(gene_id_GRCh38, pTriplo_liftover_GRCh37))
names(duplicate_dosage) <- c("GeneID", "pTriplo")
combined <- dplyr::inner_join(duplicate_dists, duplicate_dosage, by = "GeneID")
combined$pTriploAdjust <- pmin(pmax(combined$pTriplo,0.001), 0.999)
combined$pTriploAdjust <- log(combined$pTriploAdjust / (1 - combined$pTriploAdjust))
combined$TE_lengthAdjust <- log(combined$TE_length)
combined$pTriploHigh <- as.integer(combined$pTriplo >= 0.5)

ggplot(combined, aes(x = pTriploAdjust)) +
  geom_histogram()

ggplot(combined, aes(x=pTriplo)) +
  geom_histogram()

ggplot(combined, aes(x = TE_length)) +
  geom_histogram()

ggplot(combined, aes(x=TE_lengthAdjust)) +
  geom_histogram()

ml <- zeroinfl(Distance ~ pTriploHigh, data = combined, dist = "negbin")
summary(ml)

ml <- lm(Distance ~ pTriploAdjust, data = combined)
summary(ml)
