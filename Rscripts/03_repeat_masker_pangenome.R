library(dplyr)
library(readr)
library(stringr)
library(scales)
library(tidyr)
library(ggplot2)

setwd("/Users/polenpineda/Documents/buffalo_pangenome/2025/final_analysis/vcfwave_variations/repeat_masker/")
system("tail -n +3 clean_SV_50bp_filtered.fa.out | grep -v '^$' | sed -e 's/  */\\t/g' | sed 's/^\\t//' > clean_SV_50bp_filtered.txt")

repm <- "clean_SV_50bp_filtered.txt"
repm <- read_tsv(repm, col_names = c("score", "divergence", "deletion", "insertion", "query_name", "query_begin", 
                                        "query_end", "query_left", "strand","repeat", "family", "repeat_begin", 
                                        "repeat_end","repeat_left","ID", "remarks"))
repm_df <- repm %>% separate(query_name, 
                                into = c("chr", "seq_start", "copy_ID", "allele","type"), 
                                sep = "-", 
                                remove = FALSE)

repm_df_filt <- repm_df %>% 
  mutate(query_align_len = query_end - query_begin + 1, perc_id = 1 - (divergence/100)) %>%
  filter(perc_id > 0.6) %>%
  mutate(class = family, size = query_align_len) %>%
  separate(class, into = c("class","subclass"), sep = "/")

repm_df_subset <- repm_df_filt %>%
  select(ID, type, chr, seq_start, size, allele, `repeat`, family, query_align_len) %>%
  mutate(seq_start = as.numeric(seq_start),
         size = as.numeric(size),
         seq_end = seq_start + size)

repm_df_subset$size_bin <- cut(
  repm_df_subset$size,
  breaks = c(50, 1000, 100000, 500000, Inf),
  labels = c("50-1,000", "1,000-100,000", "100,000-500,000", ">500,000"),
  include.lowest = TRUE)

## bubble plot of the indel sizes
plot_chr_indel <- repm_df_subset %>%
  ggplot(aes(x = seq_start/1e6, y = factor(chr, levels = as.character(1:23)), size = size_bin, color = type)) +
  geom_point(alpha = 0.3) +
  scale_size_manual(
    values = c("50-1,000" = 0.5, "1,000-100,000" = 1, "100,000-500,000" = 3, ">500,000" = 6),
    name = "Size Range") +
  scale_x_continuous(
    breaks = seq(0, 300, by = 50),
    labels = scales::comma_format(),
    limits = c(0, 300)) +
  labs(x = "Genomic Position (Mb)", y = "Chromosome", title = "Structural variations (InDel)") +
  theme_classic()
plot_chr_indel
ggsave(plot_chr_indel, file = "plot_chr_indel.tiff", width = 8, height = 4, dpi = 600)

## separate bubble plot del and ins

# Bubble plot for Deletions (DEL)
plot_chr_del <- repm_df_subset %>%
  filter(type == "del") %>%
  ggplot(aes(x = seq_start/1e6, y = factor(chr, levels = as.character(1:23)), size = size_bin)) +
  geom_point(alpha = 0.3, color = "red") +  # Color for deletions
  scale_size_manual(
    values = c("50-1,000" = 0.5, "1,000-100,000" = 1, "100,000-500,000" = 3, ">500,000" = 6),
    name = "Size Range"
  ) +
  scale_x_continuous(
    breaks = seq(0, 300, by = 50),
    labels = scales::comma_format(),
    limits = c(0, 300)
  ) +
  labs(x = "Genomic Position (Mb)", y = "Chromosome", title = "Structural Variations (Deletions)") +
  theme_classic()

# Bubble plot for Insertions (INS)
plot_chr_ins <- repm_df_subset %>%
  filter(type == "ins") %>%
  ggplot(aes(x = seq_start/1e6, y = factor(chr, levels = as.character(1:23)), size = size_bin)) +
  geom_point(alpha = 0.3, color = "blue") +  # Color for insertions
  scale_size_manual(
    values = c("50-1,000" = 0.5, "1,000-100,000" = 1, "100,000-500,000" = 3, ">500,000" = 6),
    name = "Size Range"
  ) +
  scale_x_continuous(
    breaks = seq(0, 300, by = 50),
    labels = scales::comma_format(),
    limits = c(0, 300)
  ) +
  labs(x = "Genomic Position (Mb)", y = "Chromosome", title = "Structural Variations (Insertions)") +
  theme_classic()

# Display the plots
plot_chr_del
plot_chr_ins

summary_repm_family <- repm_df_subset %>%
  group_by(family) %>%
  summarise(count = n(), size = sum(query_align_len))

summary_type <- repm_df_filt %>%
  group_by(class) %>%
  summarise(bp = sum(query_align_len))

SV_type <- repm_df_filt %>%
  group_by(type) %>%
  summarise(bp = sum(query_align_len), count = n())

rna_types <- c("snRNA", "rRNA", "tRNA", "srpRNA", "scRNA")
summary_type_edited <- repm_df_filt %>%
  mutate(class = ifelse(class %in% rna_types, "RNA", class)) %>%
  group_by(class) %>%
  summarise(bp = sum(query_align_len)) %>%
  ungroup()
  
plot_summary_type <- summary_type_edited %>%
  ggplot(aes(x = reorder(class, -bp), y = bp/1e6, fill = class)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  scale_y_continuous(labels = comma) +
  labs(y = "Size (Mbp)", x = "") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
plot_summary_type

ggsave(plot_summary_type, file = "plot_summary_type.png", width = 3.5, height = 4, dpi = 300)

#### filter for the highest peak

repm_df_filt_small <- repm_df_filt %>%
  filter(query_align_len > 100 & query_align_len < 250) %>%
  group_by(query_align_len, `repeat`) %>%
  summarise(count = n())

repm_df_filt_big <- repm_df_filt %>%
  filter(query_align_len > 1100 & query_align_len < 1400) %>%
  group_by(query_align_len, `repeat`) %>%
  summarise(count = n())

