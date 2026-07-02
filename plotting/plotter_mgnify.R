library(dplyr)
library(patchwork)
library(purrr)
library(stringr)
library(ggplot2)
library(scales)

setwd("//wsl.localhost/Ubuntu/home/adria/tfm_app")
df_mgnify <- read.table("MGNIFY_FINAL_GAD.tsv", header=TRUE, sep="\t")

completes <- df_mgnify %>% 
  filter(Status=="COMPLETE") %>% 
  filter(Domain=="Bacteria")

nobact <- df_mgnify %>% 
  filter(Status=="COMPLETE") %>% 
  filter(Domain!="Bacteria")

biome_order <- df_mgnify %>% 
  group_by(Biome) %>% 
  summarise(total = n()) %>% 
  arrange(total) %>% 
  pull(Biome)

df_mgnify$Biome <- factor(df_mgnify$Biome, levels = biome_order)

proporciones <- ggplot(df_mgnify, aes(x = Biome, fill = Status)) +
  geom_bar(position = "fill") + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) + 
  scale_fill_manual(values = c("COMPLETE" = "lightgreen", 
                               "PARTIAL_ENZYME_ONLY" = "lightblue", 
                               "PARTIAL_TRANSPORTER_ONLY" = "lightyellow", 
                               "ABSENT" = "salmon")) +
  coord_flip() + 
  labs(title = "GAD system distribution in MGnify",
       subtitle = "Relative proportion by Biome",
       x = "Biome", y = "Proportion (%)", fill = "GAD status") +
  theme_minimal()

absolutos <- ggplot(df_mgnify, aes(x = Biome)) +
  geom_bar(fill = "grey80") +
  geom_text(stat = "count", aes(label = after_stat(count)), hjust = -0.1, size = 3) +
  coord_flip() +
  labs(title = "Genomic representation",
       subtitle = paste("Total MAGs:", nrow(df_mgnify)),
       x = "", y = "Total MAGs") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )

panel_mgnify <- proporciones + absolutos + 
  plot_layout(widths = c(2, 1), guides = "collect") & 
  theme(legend.position = "bottom")

panel_mgnify
ggsave("results/mgnify_biome_distribution.png", panel_mgnify, width = 20, height = 12)


library(tidyverse)
library(patchwork)

df_mgnify_unified <- df_mgnify %>%
  mutate(
    Status = case_when(
      Status == "PARTIAL_ENZYME_ONLY" ~ "PARTIAL_ENZYME",
      Status == "PARTIAL_TRANSPORTER_ONLY" ~ "PARTIAL_TRANSPORTER",
      TRUE ~ Status
    ),
    Phylum_Unified = case_when(
      str_detect(Phylum, "Firmicutes|Bacillota") ~ "Bacillota / Firmicutes",
      str_detect(Phylum, "Proteobacteria|Pseudomonadota") ~ "Pseudomonadota / Proteobacteria",
      str_detect(Phylum, "Actinobacteriota|Actinomycetota") ~ "Actinomycetota",
      str_detect(Phylum, "Bacteroidetes|Bacteroidota") ~ "Bacteroidota",
      is.na(Phylum) | Phylum == "" | nchar(as.character(Phylum)) <= 1 ~ "Unknown",
      TRUE ~ as.character(Phylum)
    )
  )

main_phyla <- df_mgnify_unified %>%
  count(Phylum_Unified) %>%
  mutate(perc = n / sum(n) * 100) %>%
  filter(perc >= 1.0 & Phylum_Unified != "Unknown") %>% 
  pull(Phylum_Unified)

df_mgnify_plot <- df_mgnify_unified %>%
  mutate(Phylum_Clean = ifelse(Phylum_Unified %in% main_phyla, Phylum_Unified, "Others")) %>%
  filter(Phylum_Clean != "Unknown")

phylum_order <- df_mgnify_plot %>%
  count(Phylum_Clean) %>%
  arrange(n) %>%
  pull(Phylum_Clean)

df_mgnify_plot$Phylum_Clean <- factor(df_mgnify_plot$Phylum_Clean, levels = phylum_order)

p_prop <- ggplot(df_mgnify_plot, aes(x = Phylum_Clean, fill = Status)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_fill_manual(values = c("COMPLETE" = "lightgreen", 
                               "PARTIAL_ENZYME" = "lightblue", 
                               "PARTIAL_TRANSPORTER" = "lightyellow", 
                               "ABSENT" = "salmon")) +
  coord_flip() +
  labs(title = "GAD System Distribution (MGnify)", x = "Phylum", y = "Relative Proportion") +
  theme_minimal()

p_abs <- ggplot(df_mgnify_plot, aes(x = Phylum_Clean)) +
  geom_bar(fill = "grey80") +
  geom_text(stat = "count", aes(label = after_stat(count)), hjust = -0.1, size = 3) +
  coord_flip() +
  labs(title = "Total MAGs", x = "", y = "Count") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

panel_distribucion <- p_prop + p_abs + 
  plot_layout(widths = c(2, 1), guides = "collect") & 
  theme(legend.position = "bottom")

ggsave("results/MGnify_Phylum_Distribution.png", panel_distribucion, width = 24, height = 8, bg = "white")

enrich_data <- df_mgnify_plot %>%
  group_by(Phylum_Clean) %>%
  summarise(COMPLETE = sum(Status == "COMPLETE"), Total = n(), .groups = 'drop') %>%
  rename(Phylum = Phylum_Clean)

t_comp_g <- sum(enrich_data$COMPLETE)
t_gen_g  <- sum(enrich_data$Total)

enrich_res <- enrich_data %>%
  rowwise() %>%
  mutate(
    a = COMPLETE, b = t_comp_g - COMPLETE,
    c = Total - COMPLETE, d = (t_gen_g - t_comp_g) - (Total - COMPLETE),
    fisher_res = list(fisher.test(matrix(c(a, c, b, d), nrow = 2))),
    p_val = fisher_res$p.value,
    OR = as.numeric(fisher_res$estimate),
    log2OR = ifelse(OR == 0, log2(0.01), log2(OR))
  ) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p_val, method = "fdr")) %>%
  arrange(log2OR) %>%
  mutate(Phylum = factor(Phylum, levels = Phylum))

plot_enrich_phylum <- ggplot(enrich_res, aes(x = Phylum, y = log2OR, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  geom_text(aes(label = ifelse(fdr < 0.05, "*", "")), 
            hjust = ifelse(enrich_res$log2OR > 0, -0.5, 1.5), size = 8, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(nrow(enrich_res))) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Phylum Enrichment in MGnify (Complete Systems)",
       subtitle = "Log2 Odds Ratio | Asterisk (*): FDR < 0.05",
       x = "Phylum", y = "Log2 Odds Ratio", fill = "Phylum") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_text(face = "bold", size = 10))

ggsave("results/MGnify_Phylum_Enrichment.png", plot_enrich_phylum, width = 12, height = 8, bg = "white")




####### BIOMAS ###########

library(tidyverse)

biomes_to_keep <- df_mgnify %>%
  count(Biome) %>%
  filter(n >= 10) %>%
  pull(Biome)

df_biome_enrich_prep <- df_mgnify %>%
  mutate(Status = case_when(
    Status == "PARTIAL_ENZYME_ONLY" ~ "PARTIAL_ENZYME",
    Status == "PARTIAL_TRANSPORTER_ONLY" ~ "PARTIAL_TRANSPORTER",
    TRUE ~ Status
  )) %>%
  filter(Biome %in% biomes_to_keep) %>%
  group_by(Biome) %>%
  summarise(
    COMPLETE = sum(Status == "COMPLETE"),
    Total = n(),
    .groups = 'drop'
  )

total_comp_mg <- sum(df_biome_enrich_prep$COMPLETE)
total_gen_mg <- sum(df_biome_enrich_prep$Total)

enrich_biome_res <- df_biome_enrich_prep %>%
  rowwise() %>%
  mutate(
    a = COMPLETE,
    b = total_comp_mg - COMPLETE,
    c = Total - COMPLETE,
    d = (total_gen_mg - total_comp_mg) - (Total - COMPLETE),
    
    fisher_res = list(fisher.test(matrix(c(a, c, b, d), nrow = 2))),
    p_val = fisher_res$p.value,
    OR = as.numeric(fisher_res$estimate),
    log2OR = ifelse(OR == 0, log2(0.01), log2(OR))
  ) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p_val, method = "fdr")) %>%
  arrange(log2OR) %>%
  mutate(Biome = factor(Biome, levels = Biome))

plot_biome_enrich <- ggplot(enrich_biome_res, aes(x = Biome, y = log2OR, fill = Biome)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  geom_text(aes(label = ifelse(fdr < 0.05, "*", "")), 
            hjust = ifelse(enrich_biome_res$log2OR > 0, -0.5, 1.5), 
            size = 7, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(nrow(enrich_biome_res))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  labs(
    title = "Biome Enrichment in MGnify",
    subtitle = "Fisher's Exact Test (Complete Systems) | Asterisk (*): FDR < 0.05",
    x = "Biome", 
    y = "Log2 Odds Ratio (Enrichment)",
    fill = "Biome Type"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 7)
  )

ggsave("results/MGnify_Biome_Enrichment_Biome.png", plot_biome_enrich, width = 12, height = 10, bg = "white")

print(plot_biome_enrich)

