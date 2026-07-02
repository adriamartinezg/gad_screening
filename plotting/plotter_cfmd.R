library(dplyr)
library(patchwork)
library(ggplot2)
library(forcats)

setwd("//wsl.localhost/Ubuntu/home/adria/tfm_app")
df_raw <- read.delim("GAD_SYSTEM_MAG_LEVEL_RESULTS.tsv", sep="\t")

df <- df_raw %>%
  filter(category != "" & !is.na(category)) %>%
  filter(mapply(grepl, dataset_name, MAG_id))
category_order <- df %>% 
  group_by(category) %>% 
  summarise(total = n()) %>% 
  arrange(total) %>% 
  pull(category)

df$category <- factor(df$category, levels = category_order)

proporciones <- ggplot(df, aes(x = category, fill = status)) +
  geom_bar(position = "fill") + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) + 
  scale_fill_manual(values = c("COMPLETE" = "lightgreen", 
                               "PARTIAL_ENZYME" = "lightblue", 
                               "PARTIAL_TRANSPORTER" = "lightyellow", 
                               "ABSENT" = "salmon")) +
  coord_flip() + 
  labs(title = "GAD system distribution",
       subtitle = "Relative proportion by category",
       x = "Category", y = "Proportion (%)", fill = "GAD status") +
  theme_minimal()

absolutos <- ggplot(df, aes(x = category)) +
  geom_bar(fill = "grey80") +
  geom_text(stat = "count", aes(label = after_stat(count)), hjust = -0.1, size = 3) +
  coord_flip() +
  labs(title = "Genomic representation",
       subtitle = paste("Total samples:", nrow(df)),
       x = "", y = "Total samples") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )

panel_final <- proporciones + absolutos + 
  plot_layout(widths = c(2, 1), guides = "collect") & 
  theme(legend.position = "bottom")

panel_final

ggsave("results/cfmd_proportions_categories.png", panel_final, width = 25, height = 10)


df_tax <- read.delim("final_taxonomy_results_complete.tsv", sep = "\t")

df_tax_clean <- df_tax %>%
  filter(System_Status %in% c("COMPLETE", "ABSENT", "PARTIAL_ENZYME", "PARTIAL_TRANSPORTER")) %>%
  filter(!is.na(Phylum) & Phylum != "" & nchar(as.character(Phylum)) > 1) %>%
  filter(Kuperkingdom=="Bacteria")
  droplevels()

dfc <- df_tax_clean %>%
  group_by(Phylum) %>%
  summarise(
    COMPLETE = sum(System_Status == "COMPLETE"),
    Total_Genomes = n(),
    .groups = 'drop'
  )

total_complete_global <- sum(dfc$COMPLETE)
total_genomes_global <- sum(dfc$Total_Genomes)

enrich_df <- dfc %>%
  rowwise() %>%
  mutate(
    a = COMPLETE,
    b = total_complete_global - COMPLETE,
    c = Total_Genomes - COMPLETE,
    d = (total_genomes_global - total_complete_global) - (Total_Genomes - COMPLETE),
    
    fisher_result = list(fisher.test(matrix(c(a, c, b, d), nrow = 2))),
    p_val = fisher_result$p.value,
    OR = as.numeric(fisher_result$estimate),
    log2OR = ifelse(OR == 0, log2(0.01), log2(OR)) 
  ) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p_val, method = "fdr"))


enrich_plot_data <- enrich_df %>%
  arrange(log2OR) %>%
  mutate(Phylum = factor(Phylum, levels = Phylum))

phylum_order_abs <- df_tax_clean %>% 
  count(Phylum) %>% 
  arrange(n) %>% 
  pull(Phylum)

df_tax_clean$Phylum <- factor(df_tax_clean$Phylum, levels = phylum_order_abs)

p_prop <- ggplot(df_tax_clean, aes(x = Phylum, fill = System_Status)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_fill_manual(values = c("COMPLETE" = "lightgreen", 
                               "PARTIAL_ENZYME" = "lightblue", 
                               "PARTIAL_TRANSPORTER" = "lightyellow", 
                               "ABSENT" = "salmon")) +
  coord_flip() +
  labs(title = "GAD System Distribution (cFMD)", x = "Phylum", y = "Relative Proportion") +
  theme_minimal()

p_abs <- ggplot(df_tax_clean, aes(x = Phylum)) +
  geom_bar(fill = "grey80") +
  geom_text(stat = "count", aes(label = after_stat(count)), hjust = -0.1, size = 3.5) +
  coord_flip() +
  labs(title = "Total MAGs per Phylum", x = "", y = "Count") +
  theme_classic() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank())

panel_dist_phylum <- p_prop + p_abs + 
  plot_layout(widths = c(2, 1), guides = "collect") & 
  theme(legend.position = "bottom")

ggsave("results/cFMD_Phylum_Distribution_Panel.png", panel_dist_phylum, width = 25, height = 8, bg = "white")

p_enrich_solo <- ggplot(enrich_plot_data, aes(x = Phylum, y = log2OR, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  geom_text(aes(label = ifelse(fdr < 0.05, "*", "")), 
            hjust = ifelse(enrich_plot_data$log2OR > 0, -0.5, 1.5), 
            size = 8, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(nrow(enrich_plot_data))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  labs(title = "Phylum Enrichment in cFMD (Complete Systems)",
       subtitle = "Log2 Odds Ratio | Asterisk (*): FDR < 0.05",
       x = "Phylum", y = "Log2 Odds Ratio", fill = "Phylum") +
  theme_minimal() +
  theme(
    legend.position = "right", 
    axis.text.y = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank()
  )

ggsave("results/cFMD_Phylum_Enrichment.png", p_enrich_solo, width = 12, height = 8, bg = "white")





df_results <- read.delim("GAD_SYSTEM_MAG_LEVEL_RESULTS.tsv", sep="\t")

df_cat_enrich <- df_results %>%
  dplyr::filter(status %in% c("COMPLETE", "ABSENT", "PARTIAL_ENZYME", "PARTIAL_TRANSPORTER")) %>%
  dplyr::filter(!is.na(category) & category != "") %>%
  group_by(category) %>%
  summarise(COMPLETE = sum(status == "COMPLETE"),
            Total = n(), .groups = 'drop')

t_comp_cat_g <- sum(df_cat_enrich$COMPLETE)
t_gen_cat_g <- sum(df_cat_enrich$Total)

enrich_cat_plot_data <- df_cat_enrich %>%
  rowwise() %>%
  mutate(
    a = COMPLETE,
    b = t_comp_cat_g - COMPLETE,
    c = Total - COMPLETE,
    d = (t_gen_cat_g - t_comp_cat_g) - (Total - COMPLETE),
    fisher_res = list(fisher.test(matrix(c(a, c, b, d), nrow = 2))),
    p_val = fisher_res$p.value,
    OR = as.numeric(fisher_res$estimate),
    log2OR = ifelse(OR == 0, log2(0.01), log2(OR))
  ) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p_val, method = "fdr")) %>%
  arrange(log2OR) %>%
  mutate(category = factor(category, levels = category))

plot_cat_solo <- ggplot(enrich_cat_plot_data, aes(x = category, y = log2OR, fill = category)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  geom_text(aes(label = ifelse(fdr < 0.05, "*", "")), 
            hjust = ifelse(enrich_cat_plot_data$log2OR > 0, -0.5, 1.5), 
            size = 7, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(nrow(enrich_cat_plot_data))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  labs(title = "Enrichment of GAD System by Environment", 
       subtitle = "Log2 Odds Ratio (Complete Systems) | Asterisk (*): FDR < 0.05",
       x = "Category", 
       y = "Log2 Odds Ratio (Enrichment > 0 / Depletion < 0)",
       fill = "Environment Category") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"), 
    legend.position = "right",                           
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank()
  )

ggsave("results/cfmd_Category_Enrichment.png", plot_cat_solo, width = 12, height = 8, bg = "white")

print(plot_cat_solo)