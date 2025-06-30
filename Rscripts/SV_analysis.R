library(readr)
library(dplyr)
library(ggplot2)
library(scales)
library(stringr)
library(ggtext)
library(UpSetR)

SV <- "/Users/polenpineda/Documents/buffalo_pangenome/2025/final_analysis/SV_analysis/SV_50bp_filtered.vcf.length.vcf"
SV_df <- read_delim(SV)

SV_vcf <- "/Users/polenpineda/Documents/buffalo_pangenome/2025/final_analysis/SV_analysis/SV_50bp_filtered.vcf"
SV_vcf <- read_delim(SV_vcf, skip = 50)

SV_vep <- "/Users/polenpineda/Documents/buffalo_pangenome/2025/final_analysis/vcfwave_variations/variant_effect/SV/SV_50bp_filtered.vcf.res.txt"
SV_vep <- read_delim(SV_vep, skip = 30)

### Variant Effect Predictor (VEP)

SV_vep <- SV_vep %>% mutate(ID = `#Uploaded_variation`)
SV_vep_unique <- SV_vep %>%
  mutate(priority = case_when(
    grepl("HIGH", Extra) ~ 1,
    grepl("MODERATE", Extra) ~ 2,
    grepl("LOW", Extra) ~ 3,
    grepl("MODIFIER", Extra) ~ 4,
    TRUE ~ 5
  )) %>%
  arrange(ID, priority) %>%
  distinct(ID, .keep_all = TRUE) %>%
  select(-priority)

SV_vep_vcf <- merge(SV_vcf, SV_vep_unique, by = "ID", all.x = TRUE)
SV_vep_vcf <- SV_vep_vcf %>%
  mutate(impact = str_extract(Extra, "(?<=IMPACT=)[^;]+"))

summarise_impact <- SV_vep_vcf %>% group_by(impact) %>% summarise(count = n())

### Impact and size
## bubble plot
### SV_sizes

SV_vep_vcf$ref_size <- nchar(SV_vep_vcf$REF)
SV_vep_vcf$alt_size <- nchar(SV_vep_vcf$ALT)
SV_vep_vcf$sv_size <- abs(SV_vep_vcf$ref_size - SV_vep_vcf$alt_size)
#SV_vep_vcf$sv_size <- abs(SV_vep_vcf$REF - SV_vep_vcf$ALT)


SV_indel <- SV_vep_vcf %>% filter(ref_size == 1 | alt_size == 1)
SV_others <- SV_vep_vcf %>% filter(ref_size > 1 & alt_size > 1)

sum_all_SV <- sum(SV_vep_vcf$sv_size) ## repeat masker 278736649 bp
sum_all_SV
  
SV_indel$size_bin <- cut(
  SV_indel$sv_size,
  breaks = c(0, 50, 1000, 10000, 100000, Inf),
  labels = c("<50bp", "50bp-1kb", "1kb-10kb", "10kb-100kb", ">100kb"),
  right = FALSE
)

SV_others$size_bin <- cut(
  SV_others$sv_size,
  breaks = c(0, 50, 1000, 10000, 100000, Inf),
  labels = c("<50bp", "50bp-1kb", "1kb-10kb", "10kb-100kb", ">100kb"),
  right = FALSE
)

sample_cols <- grep("^SW|^RV", names(SV_indel), value = TRUE)

SV_indel_AF <- SV_indel %>%
  rowwise() %>%
  mutate(
    alt_count = sum(c_across(all_of(sample_cols)) == "1", na.rm = TRUE),
    called_genotypes = sum(c_across(all_of(sample_cols)) %in% c("0", "1")),
    allele_freq = alt_count / (called_genotypes * 2)) %>%
  ungroup() %>%
  mutate(af_bin = cut(
      allele_freq,
      breaks = c(0, 0.05, 0.1, 0.5, 1),
      labels = c("<0.05", "0.05-0.1", "0.1-0.5", ">0.5"),
      include.lowest = TRUE,
      right = FALSE))

SV_others_AF <- SV_others %>%
  rowwise() %>%
  mutate(
    alt_count = sum(c_across(all_of(sample_cols)) == "1", na.rm = TRUE),
    called_genotypes = sum(c_across(all_of(sample_cols)) %in% c("0", "1")),
    allele_freq = alt_count / (called_genotypes * 2)) %>%
  ungroup() %>%
  mutate(af_bin = cut(
    allele_freq,
    breaks = c(0, 0.05, 0.1, 0.5, 1),
    labels = c("<0.05", "0.05-0.1", "0.1-0.5", ">0.5"),
    include.lowest = TRUE,
    right = FALSE))

SV_ins_AF <- SV_indel_AF %>%
  filter(ref_size==1)

SV_del_AF <- SV_indel_AF %>%
  filter(alt_size==1)

sum(SV_ins_AF$sv_size) 
sum(SV_del_AF$sv_size)
sum(SV_others$sv_size)

# SV_ins_bubble <- SV_ins_AF %>%
#   mutate(allele_freq = as.numeric(allele_freq)) %>%
#   group_by(size_bin, af_bin) %>%
#   summarise(variant_count = n(), .groups = "drop")
SV_ins_bubble <- SV_ins_AF %>%
  mutate(
    allele_freq = as.numeric(allele_freq),
    is_high = ifelse(impact == "HIGH", 1, 0)
  ) %>%
  group_by(size_bin, af_bin) %>%
  summarise(
    variant_count = n(),
    high_count = sum(is_high),
    high_ratio = high_count / variant_count,
    .groups = "drop"
  )

SV_del_bubble <- SV_del_AF %>%
  mutate(
    allele_freq = as.numeric(allele_freq),
    is_high = ifelse(impact == "HIGH", 1, 0)
  ) %>%
  group_by(size_bin, af_bin) %>%
  summarise(
    variant_count = n(),
    high_count = sum(is_high),
    high_ratio = high_count / variant_count,
    .groups = "drop"
  )

SV_others_bubble <- SV_others_AF %>%
  mutate(
    allele_freq = as.numeric(allele_freq),
    is_high = ifelse(impact == "HIGH", 1, 0)
  ) %>%
  group_by(size_bin, af_bin) %>%
  summarise(
    variant_count = n(),
    high_count = sum(is_high),
    high_ratio = high_count / variant_count,
    .groups = "drop"
  )

## combined

# Add a new column to indicate SV type
SV_ins_bubble <- SV_ins_bubble %>% mutate(sv_type = "Insertion")
SV_del_bubble <- SV_del_bubble %>% mutate(sv_type = "Deletion")
SV_others_bubble <- SV_others_bubble %>% mutate(sv_type = "Others")

# Combine the two datasets
# SV_bubble_combined <- bind_rows(SV_ins_bubble, SV_del_bubble, SV_others_bubble)
SV_bubble_combined <- bind_rows(SV_ins_bubble, SV_del_bubble)

# Plot with facet_wrap
plot_SV_bubble <- ggplot(SV_bubble_combined, aes(x = af_bin, y = size_bin, size = variant_count, fill = high_ratio)) +
  geom_point(alpha = 0.7, shape = 21, color = "white") +
  geom_text(aes(label = comma(variant_count)), color = "black", size = 3, vjust = -2) +  # label above point
  facet_wrap(~ sv_type) +
  scale_size_area(max_size = 12, labels = comma) +
  scale_fill_gradient(low = "lightgreen", high = "red", name = "High Impact Ratio") +
  theme_classic() +
  labs(
    x = "Allele Frequency",
    y = "SV Sizes",
    size = "Variant Count",
  )
plot_SV_bubble
ggsave(plot_SV_bubble, file = "/Users/polenpineda/Documents/buffalo_pangenome/2025/final_analysis/SV_analysis/plot_SV_bubble.png",
       width = 6, height = 4, dpi = 300)

## more details
plot_SV_bubble_detailed <- ggplot(SV_bubble_combined, aes(x = af_bin, y = size_bin, size = variant_count, fill = high_ratio)) +
  geom_point(alpha = 0.7, shape = 21, color = "black") +
  ggtext::geom_richtext(
    aes(label = paste0("<span style='color:red;'>", high_count, "</span>/<span style='color:black;'>", comma(variant_count), "</span>")),
    fill = NA, label.color = NA, size = 3, vjust = -1
  ) +
  facet_wrap(~ sv_type) +
  scale_size_area(max_size = 12, labels = comma) +
  scale_fill_gradient(low = "lightgreen", high = "tomato", name = "High Impact Ratio") +
  theme_classic() +
  labs(
    x = "Allele Frequency",
    y = "SV Size Bin",
    size = "Variant Count",
    title = "SV Allele Frequency vs. Size (Bubble Plot)"
  )
plot_SV_bubble_detailed

ggsave(plot_SV_bubble_detailed, file = "/Users/polenpineda/Documents/buffalo_pangenome/2025/final_analysis/SV_analysis/plot_SV_bubble_detailed.png",
       width = 8, height = 5, dpi = 300)


### genes

gtf <- read_delim("/Users/polenpineda/Documents/buffalo_pangenome/2025/variant_effect_predictor/PCC_UOA_SB_1v2_genomic.gtf",
                  comment = "#", col_names = FALSE)
gtf_geneid_transcript <- gtf %>%
  mutate(
    transcript_id = gsub('"', '', str_match(X9, "transcript_id\\s+([^;]+)")[,2]),
    Gene = gsub('"', '', str_match(X9, "GeneID\\:+([^;]+)")[,2]),
    protein_id = gsub('"', '', str_match(X9, "protein_id\\s+([^;]+)")[,2]),
    gene_name = gsub('"', '', str_match(X9, "gene_id\\s+([^;]+)")[,2])
  ) %>%
  select(Gene, gene_name) %>%
  distinct() %>%
  na.omit()

SV_vep_vcf_gene <- merge(SV_vep_vcf, gtf_geneid_transcript, by = "Gene", all.x = TRUE)

### what genes are affected by the high impact variants?

## high impact variants >100Kb
SV_high_impact_100kb <- SV_vep_vcf_gene %>% filter(impact == "HIGH", sv_size > 100000) %>%
  select(-REF, -ALT, -INFO, -Allele)

SV_high_impact <- SV_vep_vcf_gene %>% filter(impact == "HIGH")

### which were distinct in swamp and river?
sw_cols <- c("SWCU#0", "SWWA#0")
rv_cols <- grep("^RV", names(SV_high_impact), value = TRUE)

swamp_river_distinct <- SV_high_impact[
  (SV_high_impact$`SWCU#0` == "0" & SV_high_impact$`SWWA#0` == "0") & 
    apply(SV_high_impact[ , rv_cols], 1, function(x) all(x %in% c("1", "."))),]

### get distinct for each assemblies

# Combine all assemblies (swamp and river)
all_assemblies <- c(sw_cols, rv_cols)

for (assembly in all_assemblies) {
  # Exclude the current assembly from others
  other_assemblies <- setdiff(all_assemblies, assembly)
  
  # Filter for distinct variants
  distinct_variants <- SV_high_impact[
    SV_high_impact[[assembly]] %in% c("1", ".") &
      apply(SV_high_impact[, other_assemblies], 1, function(x) all(x %in% c("0", "."))),
  ]
  
  cols_to_drop <- intersect(c("REF", "ALT", "Allele", "INFO"), colnames(distinct_variants))
  distinct_variants_filt <- distinct_variants %>% select(-all_of(cols_to_drop))
  
  # Create a valid variable name
  df_name <- paste0("distinct_", gsub("[^A-Za-z0-9]", "_", assembly))
  
  # Assign to environment
  assign(df_name, distinct_variants_filt)
}


## identify swamp and river genotypes
rv_cols <- grep("^RV", names(SV_vep_vcf), value = TRUE)
sw_cols <- grep("^SW", names(SV_vep_vcf), value = TRUE)

# ## upset plot
sample_cols <- c("RVAZ#1","RVAZ#2","RVCU#0","RVND#1","RVND#2","RVNR#1","RVNR#2","RVUO#0","SWCU#0","SWWA#0")

upset_gen_df <- SV_vep_vcf %>%
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


### venn diagram

library(VennDiagram)
library(grid)

# Define the areas
venn.plot <- draw.pairwise.venn(
  area1 = 8895205,            # Total graph_w_reads
  area2 = 13169118,           # Total linear_w_reads
  cross.area = 5950688,       # Shared reads
  category = c("Graph", "Linear"),
  fill = c("dodgerblue", "darkorange"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.pos = c(-20, 20),
  cat.dist = 0.05
)

grid.draw(venn.plot)

# Clear the page before plotting
grid.newpage()

# Venn diagram for second dataset
draw.pairwise.venn(
  area1 = 9482055,            # Total graph_w_reads
  area2 = 19359166,           # Total linear_w_reads
  cross.area = 5950688,       # Shared reads
  category = c("Graph", "Linear"),
  fill = c("skyblue", "tomato"),
  lty = "blank",
  cex = 2,
  cat.cex = 2,
  cat.pos = c(-20, 20),
  cat.dist = 0.05
)

library(VennDiagram)
library(gridExtra)
library(grid)
library(showtext)

# Set Arial font
font_add("Arial", "/Library/Fonts/Arial.ttf")  # Adjust path if needed
showtext_auto()

# Helper to compute and label Venn diagram with percentages
make_venn <- function(area1, area2, cross, labels, colors, title = "") {
  only1 <- area1 - cross
  only2 <- area2 - cross
  total <- only1 + only2 + cross
  
  venn.plot <- draw.pairwise.venn(
    area1 = area1,
    area2 = area2,
    cross.area = cross,
    category = labels,
    fill = colors,
    lty = "blank",
    cat.cex = 0.9,
    cex = 0.9,
    cat.fontfamily = "Arial",
    fontfamily = "Arial",
    cat.dist = 0.05,
    cat.pos = c(-20, 20),
    ind = FALSE
  )
  
  # Add % labels
  venn.plot[[7]] <- textGrob(sprintf("%.1f%%", 100 * only1 / total), gp = gpar(fontsize = 8, fontfamily = "Arial"))
  venn.plot[[8]] <- textGrob(sprintf("%.1f%%", 100 * only2 / total), gp = gpar(fontsize = 8, fontfamily = "Arial"))
  venn.plot[[9]] <- textGrob(sprintf("%.1f%%", 100 * cross / total), gp = gpar(fontsize = 8, fontfamily = "Arial"))
  
  grobTree(gList(
    textGrob(title, x = 0.5, y = 1.05, gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial")),
    venn.plot
  ))
}

# First dataset
venn1 <- make_venn(8895205, 13169118, 3796074,
                   labels = c("Graph", "Linear"),
                   colors = c("skyblue", "salmon"),
                   title = "PSB (Full Reads)")

venn2 <- make_venn(9482055, 19359166, 5950688,
                   labels = c("Graph", "Linear"),
                   colors = c("skyblue", "salmon"),
                   title = "RB (Full Reads)")

# Second dataset
venn3 <- make_venn(974911, 2049713, 308175,
                   labels = c("Graph", "Linear"),
                   colors = c("lightgreen", "gold"),
                   title = "PSB (Subset)")

venn4 <- make_venn(1077558, 2964188, 256325,
                   labels = c("Graph", "Linear"),
                   colors = c("lightgreen", "gold"),
                   title = "RB (Subset)")

# Combine all
grid.newpage()
grid.arrange(venn1, venn2, venn3, venn4,
             ncol = 2,
             top = textGrob("Graph vs Linear Read Mapping in Philippine Swamp and River Buffalo",
                            gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "Arial")))
