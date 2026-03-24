library(tidyverse)
library(scales)

source("Soil_ITS_Scripts/Jeff_proc/Soil_ITS_data_prep_includ_Jeff_up.R")


filter_reads <- read_tsv('Processed_data/CG/ITS_data_quality_summary.tsv')

filter_reads<-filter_reads %>% 
  filter(barcode %in% Samples_included)

myco_dat<-myco_dat %>% 
  filter(barcode %in% Samples_included)
Guild_dat<-Guild_dat %>% 
  filter(barcode %in% Samples_included)

# Summary stats
forward_reverse_sum <- sum(filter_reads$n_merged, na.rm = TRUE)
forward_reverse_min <- min(filter_reads$n_merged, na.rm = TRUE)
forward_reverse_max <- max(filter_reads$n_merged, na.rm = TRUE)
forward_reverse_display <- paste0(comma(forward_reverse_sum), " (", 
                                  comma(forward_reverse_min), "-", 
                                  comma(forward_reverse_max), ")")

passed_qc_sum <- sum(filter_reads$n_passed, na.rm = TRUE)
passed_qc_min <- min(filter_reads$n_passed, na.rm = TRUE)
passed_qc_max <- max(filter_reads$n_passed, na.rm = TRUE)
passed_qc_display <- paste0(comma(passed_qc_sum), " (", 
                            comma(passed_qc_min), "-", 
                            comma(passed_qc_max), ")")

assigned_fungi_sum <- sum(filter_reads$n_assigned_to_fungi_and_phylum, na.rm = TRUE)
assigned_fungi_min <- min(filter_reads$n_assigned_to_fungi_and_phylum, na.rm = TRUE)
assigned_fungi_max <- max(filter_reads$n_assigned_to_fungi_and_phylum, na.rm = TRUE)
assigned_fungi_display <- paste0(comma(assigned_fungi_sum), " (", 
                                 comma(assigned_fungi_min), "-", 
                                 comma(assigned_fungi_max), ")")

reads_after_rare_sum <- sum(filter_reads$n_after_removing_rare_otus, na.rm = TRUE)
reads_after_rare_min <- min(filter_reads$n_after_removing_rare_otus, na.rm = TRUE)
reads_after_rare_max <- max(filter_reads$n_after_removing_rare_otus, na.rm = TRUE)
reads_after_rare_display <- paste0(comma(reads_after_rare_sum), " (", 
                                   comma(reads_after_rare_min), "-", 
                                   comma(reads_after_rare_max), ")")

# OTU stats
otus_assigned_to_genus <- Guild_dat %>% 
  filter(barcode %in% Samples_included) %>% distinct(OTU) %>% nrow()
otus_assigned_to_mycorrhizal_taxa <- myco_dat %>% 
  filter(barcode %in% Samples_included) %>% distinct(OTU) %>% nrow()


# Mycorrhizal read counts per sample (display as sum/min/max)
myco_reads <- myco_dat %>% summarise(total = sum(count, na.rm = TRUE)) %>% pull(total)
perc_myco_reads <- round(100 * myco_reads / assigned_fungi_sum, 1)
perc_myco_reads_display <- paste0(perc_myco_reads, "%")

# Mycorrhizal reads per location for new row
myco_read_stats <- myco_dat%>% 
  filter(barcode %in% Samples_included) %>%
  group_by(barcode) %>%
  summarise(total_reads = sum(count, na.rm = TRUE)) %>%
  summarise(
    sum_reads = sum(total_reads, na.rm = TRUE),
    min_reads = min(total_reads, na.rm = TRUE),
    max_reads = max(total_reads, na.rm = TRUE),
    n_samples = n()
  )

myco_reads_display <- paste0(comma(myco_read_stats$sum_reads), " (", 
                             comma(myco_read_stats$min_reads), "-", 
                             comma(myco_read_stats$max_reads), ")")

# Count samples at each processing step
n_merged_samples <- sum(!is.na(filter_reads$n_merged))
n_passed_samples <- sum(!is.na(filter_reads$n_passed))
n_assigned_samples <- sum(!is.na(filter_reads$n_assigned_to_fungi_and_phylum))
n_rare_filtered_samples <- sum(!is.na(filter_reads$n_after_removing_rare_otus))
n_myco_plots <- myco_dat %>% distinct(barcode) %>% nrow()

# Final summary tibble
summary_column <- tibble(
  metric = c(
    "Merged",
    "Quality control",
    "Assigned Fungal phylum",
    "Remove rare OTUs",
    "mycorrhizal reads",
    "OTUs",
    "mycorrhizal OTUs",
    "% mycorrhizal reads"
  ),
  Soil_reads= c(
    forward_reverse_display,
    passed_qc_display,
    assigned_fungi_display,
    reads_after_rare_display,
    myco_reads_display,
    comma(otus_assigned_to_genus),
    comma(otus_assigned_to_mycorrhizal_taxa),
    perc_myco_reads_display
  ),
  samples = c(
    n_merged_samples,
    n_passed_samples,
    n_assigned_samples,
    n_rare_filtered_samples,
    myco_read_stats$n_samples,
    NA,
    NA,
    NA
  )
)

# View and write output
summary_column

write_rds(summary_column, 'Processed_data/seq_summary/soil_CG_reads_summary.rds')
