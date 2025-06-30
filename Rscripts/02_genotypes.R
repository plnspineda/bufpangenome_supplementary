library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(UpSetR)

setwd("/Users/polenpineda/Documents/buffalo_pangenome/2025/variant_effect_predictor/")
high_swamp_vep <- "/Users/polenpineda/Documents/buffalo_pangenome/2025/variant_effect_predictor/high_impact_all_chromosomes.concat.vcfwave.NORM.res.txt"
high_swamp_vep <- read_delim(high_swamp_vep, col_names = FALSE)

gtf <- read_delim("PCC_UOA_SB_1v2_genomic.gtf",
                  comment = "#", col_names = FALSE)
gtf_geneid_transcript <- gtf %>%
  mutate(
    transcript_id = gsub('"', '', str_match(X9, "transcript_id\\s+([^;]+)")[,2]),
    X4 = gsub('"', '', str_match(X9, "GeneID\\:+([^;]+)")[,2]),
    protein_id = gsub('"', '', str_match(X9, "protein_id\\s+([^;]+)")[,2]),
    gene_name = gsub('"', '', str_match(X9, "gene_id\\s+([^;]+)")[,2])
  ) %>%
  select(X4, gene_name) %>%
  distinct() %>%
  na.omit()

## create a bed file for the high_impact variant

high_df <- high_swamp_vep %>%
  separate(X2, into = c("chr","pos"), sep = ":") %>%
  separate(pos, into = c("start","end"), sep = "-")
bed <- high_df %>%
  mutate(end = if_else(is.na(end), as.numeric(start) + 1, as.numeric(end))) %>%
  select(chr, start, end) %>%
  distinct()

write_tsv(bed, file = "high_impact_variants.bed", col_names = FALSE)
id <- high_df %>% select(X1) %>% distinct()
write_tsv(id, file = "high_impact_variants.id.txt", col_names = FALSE)

### what genes are affected by the high impact variants?

high_df_gene <- merge(high_df, gtf_geneid_transcript, by = "X4", all.x = TRUE)

transcript_ablation <- high_df_gene %>%
  filter(X7 == "transcript_ablation") %>%
  select(gene_name, chr, start, end) %>%
  arrange(chr, start, end) %>%
  distinct(chr, start, end, .keep_all = TRUE)
transcript_ablation_genes <- transcript_ablation %>%
  select(gene_name) %>% distinct()

splice_acceptor_variant <- high_df_gene %>%
  filter(X7 == "splice_acceptor_variant") %>%
  select(gene_name, chr, start, end) %>%
  arrange(chr, start, end) %>%
  distinct(chr, start, end, .keep_all = TRUE)

splice_donor_variant <- high_df_gene %>%
  filter(X7 == "splice_donor_variant") %>%
  select(gene_name, chr, start, end) %>%
  arrange(chr, start, end) %>%
  distinct(chr, start, end, .keep_all = TRUE)

stop_gained <- high_df_gene %>%
  filter(X7 == "stop_gained") %>%
  select(gene_name, chr, start, end) %>%
  arrange(chr, start, end) %>%
  distinct(chr, start, end, .keep_all = TRUE)

### after extracting the high impact SNPs, find which are unique to swamp and river

vcf <- "/Users/polenpineda/Documents/buffalo_pangenome/2025/variant_effect_predictor/high_impact_variants.id.vcf"
vcf_df <- read_delim(vcf, skip = 48, col_names = TRUE)

gen_df <- vcf_df %>% select(-REF,-ALT,-QUAL,-FILTER,-INFO,-FORMAT)

## identify swamp and river genotypes
## adenyl
rv_cols <- grep("^RV", names(gen_df), value = TRUE)
sw_cols <- grep("^SW", names(gen_df), value = TRUE)

# ## upset plot
sample_cols <- c("RVAZ#1","RVAZ#2","RVCU#0","RVND#1","RVND#2","RVNR#1","RVNR#2","RVUO#0","SWCU#0","SWWA#0")

upset_gen_df <- gen_df %>%
  select(all_of(sample_cols))
#upset_gen_df[upset_gen_df == "."] <- NA
upset_gen_df[sample_cols] <- sapply(upset_gen_df[sample_cols], function(x) gsub("^\\.$", "1", x))
upset_gen_df_binary <- as.data.frame(sapply(upset_gen_df, function(x) as.numeric(x != "1")))
upset_gen_df_binary <- upset_gen_df_binary %>% mutate(`SWPC#0` = 0)

sample_cols <- c("RVAZ#1","RVAZ#2","RVCU#0","RVND#1","RVND#2","RVNR#1","RVNR#2","RVUO#0", "SWCU#0","SWWA#0", "SWPC#0")
set_colors <- ifelse(sample_cols %in% c("SWCU#0", "SWWA#0", "SWPC#0"), "orange", "navyblue")
names(set_colors) <- sample_cols


# png("plot_upset.png", width = 11, height = 5, units = "in", res = 300)
upset(upset_gen_df_binary,
      sets = sample_cols,
      order.by = "freq",
      keep.order = TRUE,
      text.scale = 0.8,
      mainbar.y.label = "Number of Variants",
      sets.x.label = "Variants per Assembly",
      sets.bar.color = set_colors)
# dev.off()

### vcf and genes

high_df_gene_subset <- high_df_gene %>% mutate(ID = X1) %>% select(ID, chr, start, end, transID = X5, effect = X7, impact = X14, gene_name)
gen_df <- gen_df %>% mutate(chr = `#CHROM`) %>% select(-`#CHROM`)
merged_info <- merge(high_df_gene_subset, gen_df, by = c("ID","chr"), all.x = TRUE)

merged_info[] <- lapply(merged_info, as.character)
sw_cols <- c("SWCU#0", "SWWA#0")
rv_cols <- grep("^RV", names(gen_df), value = TRUE)

swamp_river_distinct <- merged_info[
  (merged_info$`SWCU#0` == "0" & merged_info$`SWWA#0` == "0") & 
    apply(merged_info[ , rv_cols], 1, function(x) all(x %in% c("1", "."))),
]
swamp_river_gene_distinct <- swamp_river_distinct %>% select(gene_name) %>% distinct()
summary_swamp_river_distinct <- swamp_river_distinct %>%
  group_by(effect) %>%
  summarise(count = n())
swamp_river_gene_distinct_frameshift <- swamp_river_distinct %>%
  filter(effect == "frameshift_variant") %>% select(gene_name) %>% distinct()

spwc_distinct <- merged_info[
  (merged_info$`SWCU#0` == "1" & merged_info$`SWWA#0` == "1") & 
    apply(merged_info[ , rv_cols], 1, function(x) all(x %in% c("1", "."))),
]
spwc_distinct_gene <- spwc_distinct %>% select(gene_name) %>% distinct()

### only get SNPs

vcf_df$ref_size <- nchar(vcf_df$REF)
vcf_df$alt_size <- nchar(vcf_df$ALT)
vcf_df$sv_size <- abs(vcf_df$ref_size - vcf_df$alt_size)

SNPs_vcf <- vcf_df %>% filter(ref_size == 1 & alt_size == 1)
SNPs_indel_sm <- vcf_df %>% filter((ref_size != 1 | alt_size != 1) & abs(ref_size - alt_size) < 50)

