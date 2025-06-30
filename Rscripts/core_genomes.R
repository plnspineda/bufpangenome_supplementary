library(dplyr)
library(readr)
library(ggplot2)

setwd("/Users/polenpineda/Documents/buffalo_pangenome/2025/core_accessories_genome/")
df <- "core_genomes_stats.txt"

df <- read_delim(df)
df_core <- df %>%
  select(Chromosome, Core, Private, Shell, SWPC)

df_long_genome <- df_core %>%
  pivot_longer(cols = c(Core, Shell, Private), names_to = "Type", values_to = "Value") %>%
  mutate(Type = factor(Type, levels = c("Private", "Shell", "Core")))

df_swpc <- df_core %>%
  mutate(Type = "SWPC") %>%
  select(Chromosome, Type, Value = SWPC)

# Plot
plot_core_genomes <- ggplot() +
  geom_bar(data = df_long_genome, aes(x = factor(Chromosome), y = Value / 1e6, fill = factor(Type)),
           stat = "identity") +
  geom_bar(data = df_swpc, aes(x = factor(Chromosome), y = Value / 1e6, fill = Type),
           stat = "identity", position = position_nudge(x = 0.25), width = 0.4) +
  labs(x = "Chromosome", y = "Size (Mb)") +
  scale_fill_manual(values = c("Core" = "#427aa1", "Shell" = "#83c5be", "Private" = "#f4d35e", "SWPC" = "#d1495b")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.title = element_blank()) +
  scale_y_continuous(labels = scales::comma, limits = c(0, 350))
plot_core_genomes
ggsave(plot_core_genomes, filename = "plot_core_genomes.png", width = 8, height = 4, dpi = 300)

## make a pie chart of the variants within the pangenome graph

variant_types <- data.frame(
  type = c("SNVs", "MNPs", "indels <50bp", "indels >50bp", "SV >50bp"),
  count = c(26824142, 845153, 4674873, 172423, 7103))


variant_types$perc <- variant_types$count / sum(variant_types$count) * 100
variant_types$label <- paste0(variant_types$type, "\n", round(variant_types$per, 1), "%")
variant_types$type <- factor(variant_types$type, 
                             levels = c("SNVs", "MNPs", "indels <50bp", "indels >50bp", "SV >50bp"))

plot_variant_types <- ggplot(variant_types, aes(x = "", y = count, fill = factor(type))) +
  geom_bar(stat = "identity", width = 0.1, color = "white") +
  coord_polar("y") +
  theme_void() + 
  theme(legend.title = element_blank())
plot_variant_types
ggsave(plot_variant_types, filename = "plot_variant_types.png", width = 8, height = 4, dpi = 300)
