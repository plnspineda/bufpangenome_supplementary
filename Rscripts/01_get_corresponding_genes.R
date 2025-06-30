library(dplyr)
library(readr)
library(stringr)
library(tidyr)

setwd("/Users/polenpineda/Documents/buffalo_pangenome/2025/final_analysis/vcfwave_variations/")
getwd()

orthogroups <- read_delim("/Users/polenpineda/Documents/buffalo_pangenome/2025/variant_effect_predictor/Orthogroups.tsv")
orthogroups <- na.omit(orthogroups)

gtf <- read_delim("/Users/polenpineda/Documents/buffalo_pangenome/2025/variant_effect_predictor/PCC_UOA_SB_1v2_genomic.gtf",
                  comment = "#", col_names = FALSE)

gtf_cow <- read_delim("/Users/polenpineda/Documents/buffalo_pangenome/2025/variant_effect_predictor/GCF_002263795.3_ARS-UCD2.0_genomic.gtf",
                      comment = "#", col_names = FALSE)

SV_swamp_vep <- read_delim("/Users/polenpineda/Documents/buffalo_pangenome/2025/final_analysis/vcfwave_variations/variant_effect/SV/SV_50bp_filtered.vcf.res.txt", skip = 31, col_names = FALSE)

orthogroups_single <- orthogroups %>%
  separate(ARSUCD, into = c("cow_protein_id"), sep = ",") %>%
  separate(SWPC, into = c("protein_id"), sep = ",")

## I only have the transcipt ID in the vcf file, so I need to cross it back to get the
## corresponding geneID with the gtf file

gtf_geneid_transcript <- gtf %>%
  mutate(
    transcript_id = gsub('"', '', str_match(X9, "transcript_id\\s+([^;]+)")[,2]),
    gene_id = gsub('"', '', str_match(X9, "GeneID\\:+([^;]+)")[,2]),
    protein_id = gsub('"', '', str_match(X9, "protein_id\\s+([^;]+)")[,2])
  ) %>%
  select(transcript_id, protein_id, gene_id) %>%
  distinct() %>%
  na.omit()

geneid_unique <- unique(gtf_geneid_transcript$gene_id)

gtf_cow_geneid_transcript <- gtf_cow %>%
  mutate(
    cow_transcript_id = gsub('"', '', str_match(X9, "transcript_id\\s+([^;]+)")[,2]),
    cow_gene_id = gsub('"', '', str_match(X9, "GeneID\\:+([^;]+)")[,2]),
    cow_protein_id = gsub('"', '', str_match(X9, "protein_id\\s+([^;]+)")[,2])
  ) %>%
  select(cow_transcript_id, cow_protein_id, cow_gene_id) %>%
  distinct() %>%
  na.omit()

### getting the geneID from orthofinder
ortho_bufgene <- merge(orthogroups_single, gtf_geneid_transcript, by = "protein_id")
ortho_bufcowgene <- merge(ortho_bufgene, gtf_cow_geneid_transcript, by = "cow_protein_id")
ortho_bufcowgene <- ortho_bufcowgene %>% select(transcript_id, cow_gene_id)
ortho_dict <- merge(ortho_bufcowgene,ortho_bufgene, by = "transcript_id")

## get geneID from snpeff...
## orthofinder only have the protein id
## snpeff only have the transcript id
## i get the matched cow_gene_id and the buffalo transcript_id from orthofinder (ortho_bufcowgene)
## then i merge the swamp_snpeff transcript id to get the cow gene id

gene_dict <- merge(ortho_bufgene, ortho_bufcowgene, by = "transcript_id")

SV_swamp_vep_list <- SV_swamp_vep %>% mutate(vep = "vep") %>% select(gene_id = X4, vep)
SV_swamp_vep_dict <- merge(SV_swamp_vep_list, gene_dict, by = "gene_id")

SV_cow_equivalent_list <- SV_swamp_vep_dict %>%
  filter(gene_id != "-")

gene_list <- as.character(SV_cow_equivalent_list$cow_gene_id)
  
############

### GO ENRICHMENT
## gene enrichment
## also from this tutorial: https://learn.gencore.bio.nyu.edu/rna-seq-analysis/over-representation-analysis/
# BiocManager::install("enrichplot", force = TRUE)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("biomaRt")
#BiocManager::install(c("clusterProfiler", "org.Bt.eg.db", "AnnotationDbi"))
#BiocManager::install("AnnotationHub")
#BiocManager::install("limma")

library(AnnotationHub)
library(clusterProfiler)
library(limma)
library(org.Bt.eg.db)
library(enrichplot)

hub <- AnnotationHub()
query(hub, "org.Bt.eg.db")

swamp_go <- enrichGO(gene_list, OrgDb=org.Bt.eg.db, pvalueCutoff=1, qvalueCutoff=1, readable = TRUE)
head(swamp_go)
swamp_go_summary <- data.frame(swamp_go)
write_csv(swamp_go_summary, file = "SV_enrichment_go_ncbi.tsv", col_names = TRUE)

swamp_cluster <- groupGO(gene = gene_list, OrgDb=org.Bt.eg.db, level = 3, readable = TRUE)
swamp_cluster_summary <- data.frame(swamp_cluster) %>% filter(Count != 0)

## visualise
goplot(swamp_go)

## upset plot
upsetplot(swamp_go)

## barplot
barplot(swamp_go, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)

## dotplot
dotplot(swamp_go)

