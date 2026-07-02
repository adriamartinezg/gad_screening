library(dplyr)
library(patchwork)
library(purrr)
library(stringr)
library(ggplot2)
library(scales)
library(tidyr)

setwd("//wsl.localhost/Ubuntu/home/adria/tfm_app")
df_mgnify <- read.table("MGNIFY_FINAL_GAD.tsv", header=TRUE, sep="\t")

df_tax <- read.delim("final_taxonomy_results_complete.tsv", sep = "\t")

df_tax_clean <- df_tax %>%
  filter(System_Status %in% c("COMPLETE", "ABSENT", "PARTIAL_ENZYME", "PARTIAL_TRANSPORTER")) %>%
  filter(!is.na(Phylum) & Phylum != "" & nchar(as.character(Phylum)) > 1) %>%
  filter(Kuperkingdom=="Bacteria")

df_cfmd_prep <- df_tax_clean %>%
  filter(System_Status %in% c("COMPLETE", "ABSENT", "PARTIAL_ENZYME", "PARTIAL_TRANSPORTER")) %>%
  filter(Kuperkingdom == "Bacteria") %>%
  mutate(
    Biome = case_when(
      Category == "alcohol"                         ~ "Food: Alcohol",
      Category == "dairy"                           ~ "Food: Dairy",
      Category == "fermented_beverages"             ~ "Food: Fermented Beverages",
      Category == "fermented_fish"                  ~ "Food: Fermented Fish",
      Category == "fermented_fruits_and_vegetables" ~ "Food: Fermented Fruits/Vegetables",
      Category == "fermented_grains"                ~ "Food: Fermented Grains",
      Category == "fermented_legumes"               ~ "Food: Fermented Legumes",
      Category == "fermented_meat"                  ~ "Food: Fermented Meat",
      Category == "fermented_seeds"                 ~ "Food: Fermented Seeds",
      Category == "fermented_tubers_and_roots"       ~ "Food: Fermented Tubers/Roots",
      Category == "fish"                            ~ "Food: Fish",
      Category == "fruits_and_vegetables"           ~ "Food: Fruits/Vegetables",
      Category == "meat"                            ~ "Food: Meat",
      Category == "probiotics"                      ~ "Food: Probiotics",
      Category == "other"                           ~ "Food: Other",
      TRUE                                          ~ paste0("Food: ", str_to_title(Category))
    ),
    Status = System_Status
  ) %>%
  select(Genome_ID = MAG_id, Biome, Status)

df_mgnify_prep <- df_mgnify %>%
  filter(Domain == "Bacteria") %>% 
  select(Genome_ID, Biome, Status)

df_combined_all <- bind_rows(df_cfmd_prep, df_mgnify_prep) %>%
  filter(!is.na(Biome) & Biome != "")

df_biome_counts <- df_combined_all %>%
  group_by(Biome) %>%
  summarise(
    COMPLETE = sum(Status == "COMPLETE"),
    Total = n(),
    .groups = 'drop'
  )

t_comp_global <- sum(df_biome_counts$COMPLETE)
t_gen_global  <- sum(df_biome_counts$Total)

enrich_joint_data <- df_biome_counts %>%
  rowwise() %>%
  mutate(
    a = COMPLETE,
    b = t_comp_global - COMPLETE,
    c = Total - COMPLETE,
    d = (t_gen_global - t_comp_global) - (Total - COMPLETE),
    fisher_res = list(fisher.test(matrix(c(a, c, b, d), nrow = 2))),
    p_val = fisher_res$p.value,
    OR = as.numeric(fisher_res$estimate),
    log2OR = ifelse(OR == 0, log2(0.01), log2(OR))
  ) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p_val, method = "fdr")) %>%
  arrange(log2OR) %>%
  mutate(Biome = factor(Biome, levels = Biome))

plot_joint_final <- ggplot(enrich_joint_data, aes(x = Biome, y = log2OR, fill = Biome)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  geom_text(aes(label = ifelse(fdr < 0.05, "*", "")), 
            hjust = ifelse(enrich_joint_data$log2OR > 0, -0.5, 1.5), 
            size = 7, fontface = "bold", vjust = 0.7) +
  coord_flip() +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(nrow(enrich_joint_data))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  labs(
    title = "Enrichment of GAD System by Biome (Bacteria)", 
    subtitle = "Log2 Odds Ratio (Complete Systems vs Global Background) | Asterisk (*): FDR < 0.05",
    x = "Biome / Environment", 
    y = "Log2 Odds Ratio (Enrichment > 0 / Depletion < 0)",
    fill = "Biome Category"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9, face = "bold"), 
    legend.position = "none", 
    panel.grid.minor = element_blank()
  )

print(plot_joint_final)

ggsave("results/Enrichment_combined_biome.png", plot_joint_final, width = 12, height = 10, bg = "white")




# ==============================================================================
# INTEGRACIÓN DE REFSEQ (UTILIZANDO EL RESUMEN DE PHYLUM)
# ==============================================================================

# 1. Cargar el archivo de RefSeq por Phylum
df_refseq_raw <- read.table("RESUMEN_PHYLUM.tsv", header=TRUE, sep="\t")

df_refseq_prep <- df_refseq_raw %>%
  mutate(
    Phylum_Unified = case_when(
      str_detect(Phylum, "Firmicutes|Bacillota") ~ "Bacillota / Firmicutes",
      str_detect(Phylum, "Proteobacteria|Pseudomonadota") ~ "Pseudomonadota / Proteobacteria",
      TRUE ~ as.character(Phylum)
    ),
    Biome = "Reference (RefSeq)"
  ) %>%
  select(Biome, Phylum_Unified, COMPLETE, PARTIAL_ENZYME, PARTIAL_TRANSPORTER, ABSENT) %>%
  pivot_longer(cols = c(COMPLETE, PARTIAL_ENZYME, PARTIAL_TRANSPORTER, ABSENT), 
               names_to = "Status", values_to = "Count") %>%
  uncount(Count) %>%
  select(Biome, Phylum = Phylum_Unified, Status)

# ==============================================================================
# 2. FUSIÓN TRIPLE Y LIMPIEZA TAXONÓMICA AVANZADA
# ==============================================================================

df_triple_tax <- bind_rows(df_refseq_prep, df_cfmd_prep_tax, df_mgnify_prep_tax) %>%
  # Limpieza inicial de indeterminados
  filter(!is.na(Phylum) & Phylum != "" & Phylum != "Otros" & Phylum != "Unknown") %>%
  filter(!str_detect(Phylum, "(?i)unassigned|unclassified|unknown|assigned")) %>%
  mutate(
    Phylum_Clean = str_remove(Phylum, "_[A-Z]$")
  )

# Identificamos qué filos limpios superan el umbral de los 50 genomas totales
filos_principales <- df_triple_tax %>%
  group_by(Phylum_Clean) %>%
  tally() %>%
  filter(n >= 50) %>%
  pull(Phylum_Clean)

# Agrupamos los minoritarios bajo la etiqueta "Other Minor Phyla"
df_triple_tax_final <- df_triple_tax %>%
  mutate(
    Phylum_Final = ifelse(Phylum_Clean %in% filos_principales, Phylum_Clean, "Other Minor Phyla")
  )

# ==============================================================================
# 3. RECUENTOS Y CÁLCULO DE ENRIQUECIMIENTO GLOBAL (FISHER)
# ==============================================================================

df_enrich_tax_counts <- df_triple_tax_final %>%
  group_by(Phylum = Phylum_Final) %>%
  summarise(
    COMPLETE = sum(Status == "COMPLETE"),
    Total = n(),
    .groups = 'drop'
  )

t_comp_global_tax <- sum(df_enrich_tax_counts$COMPLETE)
t_gen_global_tax  <- sum(df_enrich_tax_counts$Total)

# test de Fisher
enrich_pure_tax_data <- df_enrich_tax_counts %>%
  rowwise() %>%
  mutate(
    a = COMPLETE,
    b = t_comp_global_tax - COMPLETE,
    c = Total - COMPLETE,
    d = (t_gen_global_tax - t_comp_global_tax) - (Total - COMPLETE),
    
    fisher_res = list(fisher.test(matrix(c(a, c, b, d), nrow = 2))),
    p_val = fisher_res$p.value,
    OR = as.numeric(fisher_res$estimate),
    log2OR = ifelse(OR == 0, log2(0.01), log2(OR))
  ) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p_val, method = "fdr"))

enrich_pure_tax_data <- enrich_pure_tax_data %>%
  arrange(Phylum == "Other Minor Phyla", log2OR) %>%
  mutate(Phylum = factor(Phylum, levels = Phylum))


# ==============================================================================
# 4. GRÁFICO DE ENRIQUECIMIENTO TAXONÓMICO 
# ==============================================================================

plot_tax_enrich_final <- ggplot(enrich_pure_tax_data, aes(x = Phylum, y = log2OR, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  geom_text(aes(label = ifelse(fdr < 0.05, "*", "")), 
            hjust = ifelse(enrich_pure_tax_data$log2OR > 0, -0.5, 1.5), 
            size = 8, fontface = "bold", vjust = 0.7) +
  coord_flip() +
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(nrow(enrich_pure_tax_data))) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  labs(
    title = "Bacterial Phylum Enrichment of the GAD System", 
    subtitle = "Global Fisher's Exact Test (RefSeq + cFMD + MGnify) | Asterisk (*): FDR < 0.05",
    x = "Bacterial Phylum", 
    y = "Log2 Odds Ratio (Enrichment > 0 / Depletion < 0)"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10, face = "bold"), 
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "none", 
    panel.grid.minor = element_blank()
  )

# Mostrar en pantalla y guardar
print(plot_tax_enrich_final)
ggsave("results/Enrichment_combined_Taxonomy_Bacteria.png", plot_tax_enrich_final, width = 11, height = 8, bg = "white")
