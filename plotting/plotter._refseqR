library(ggplot2)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(patchwork)

setwd("//wsl.localhost/Ubuntu/home/adria/tfm_app")
df<-read.table("RESUMEN_PHYLUM.tsv", header=TRUE, sep="\t")

df_long <- df %>% 
  select(Phylum, COMPLETE, PARTIAL_ENZYME, PARTIAL_TRANSPORTER, ABSENT) %>% 
  pivot_longer(cols= -Phylum, names_to ="Status", values_to = "Count")

df_long$Status <- factor(df_long$Status, 
                         levels = c("COMPLETE", "PARTIAL_ENZYME", "PARTIAL_TRANSPORTER", "ABSENT"))

proporciones <- ggplot(df_long, aes(x = reorder(Phylum, Count, sum), y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) + 
  scale_fill_manual(values = c("lightgreen", "lightblue", "lightyellow", "salmon")) +
  coord_flip() + 
  labs(title = "GAD system distribution in Bacteria",
       subtitle = "Relative proportion",
       x = "Phylum", y = "Proportion (%)", fill = "GAD state") +
  theme_minimal()

absolutos <- ggplot(df, aes(x = reorder(Phylum, Total_Genomes), y = Total_Genomes)) +
  geom_bar(stat = "identity", fill = "grey80") +
  geom_text(aes(label = Total_Genomes), hjust = -0.1, size = 3) +
  coord_flip() +
  labs(title = "Genomic representation",
       subtitle = paste("Total:", sum(df$Total_Genomes)),
       x = "", y = "Total genomes") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )

panel_final <- proporciones + absolutos + 
  plot_layout(widths = c(2, 1), guides = "collect") & 
  theme(legend.position = "bottom")

ggsave("RefSeq_plots_bact.png", panel_final, width = 20, height = 10)


dfc<-read.table("RESUMEN_CLASSES.tsv", header=TRUE, sep="\t")
dfc$perc_num <- as.numeric(gsub("%", "", dfc$"X._COMPLETE"))
dfc <- dfc %>% arrange(desc(perc_num))


dfc <- dfc %>%
  mutate(Phylum = case_when(
    Class %in% c("Bacteroidia", "Cytophagia", "Chitinophagia","Flavobacteriia", 
                 "Saprospiria", "Sphingobacteriia", "Ignavibacteria") ~ "Bacteroidota",
    Class %in% c("Bacilli", "Clostridia", "Culicoidibacteria", "Negativicutes",
                  "Erysipelotrichia", "Limnochordia", "Syntrophomonadia",
                 "Thermaerobacteria", "Thermolithobacteria") ~ "Bacillota",
    Class %in% c("Alphaproteobacteria", "Betaproteobacteria", 
                  "Gammaproteobacteria", "Deltaproteobacteria", "Epsilonproteobacteria",
                  "Zetaproteobacteria", "Hydrogenophilalia", "Oligoflexia", "Acidithiobacillales", "Magnetococcia") ~ "Pseudomonadota",
    Class %in% c("Actinomycetes", "Coriobacteriia", "Acidimicrobiia",
                  "Nitriliruptoria", "Rubrobacteria", "Thermoleophilia") ~ "Actinomycetota",
    Class %in% c("Fusobacteriia") ~ "Fusobacteriota",
    Class %in% c("Chlamydiia") ~ "Chlamydiota",
    Class %in% c("Spirochaetia") ~ "Spirochaetota",
    Class %in% c("Cyanophyceae", "Vampirovibrionophyceae") ~ "Cyanobacteriota",
    Class %in% c("Planctomycetia", "Phycisphaerae") ~ "Planctomycetota",
    Class %in% c("Anaerolineae", "Ardenticatenia", "Caldilineae", "Chloroflexia", "Dehalococcoidia",
                 "Ktedonobacteria", "Tepidiformia", "Thermoflexia") ~ "Chloroflexota",
    Class %in% c("Blastocatellia", "Holophagae","Thermoanaerobaculia", "Vicinamibacteria") ~ "Acidobacteriota",
    Class %in% c("Aquificia", "Desulfurobacteriia", "Thermosulfidibacteria") ~"Aquificota",
    Class %in% c("Atribacteria") ~ "Aquificota",
    Class %in% c("Balneolia") ~ "Balneolota",
    Class %in% c("Caldisericia") ~ "Caldisericota",
    Class %in% c("Calditrichia") ~"Calditrichota",
    Class %in% c("Chlorobia") ~ "Chlorobiota",
    Class %in% c("Chrysiogenia") ~ "Chrysiogenota",
    Class %in% c("Coprothermobacteria") ~ "Coprothermobacterota",
    Class %in% c("Deferribacteres") ~ "Deferribacterota",
    Class %in% c("Dictyoglomeria") ~ "Dyctyoglomerota",
    Class %in% c("Elusimicrobia", "Endomicrobiia") ~ "Elusimicrobiota",
    Class %in% c("Chitinispirillia", "Chitinivibrionia", "Fibrobacteria") ~ "Fibrobacterota",
    Class %in% c("Gemmatimonadetes", "Longimicrobia") ~ "Gemmatimonadota",
    Class %in% c("Kiritimatiellae", "Tichowtungiia") ~ "Kiritimatiellota",
    Class %in% c("Oligosphaeria", "Lentisphaeria") ~ "Lentisphaerota",
    Class %in% c("Bradymonadia", "Myxococcia", "Polyangiia") ~ "Myxococcota",
    Class %in% c("Nitrospinia") ~ "Nitrospinota",
    Class %in% c("Nitrospiria", "Thermodesulfovibrionia") ~ "Nitrospirota",
    Class %in% c("Rhodothermia") ~ "Rhodothermota",
    Class %in% c("Deferrisomatia", "Desulfarculia", "Desulfobaccia", "Desulfobacteria",
                 "Desulfobulbia", "Desulfomonilia", "Desulfovibrionia", "Desulfuromonadia",
                 "Dissulfuribacteria", "Syntrophia", "Syntrophobacteria", 
                 "Syntrophorhabdia", "Thermodesulfobacteria") ~ "Thermodesulfobacteriota",
    Class %in% c("Thermomicrobia") ~ "Thermomicrobiota",
    Class %in% c("Opitutia", "Verrucomicrobiia") ~ "Verrucomicrobiota",
    Class %in% c("Deinococci") ~ "Deinococcota",
    Class %in% c("Synergistia") ~ "Synergistota",
    Class %in% c("Thermotogae") ~ "Thermotogota",
    Class %in% c("Armatimonadia", "Chthonomonadia", "Fimbriimonadia") ~ "Armatimonadota",
    Class %in% c("Gracilibacteria", "Microgenomatia", "Minisyncoccia", "Patescibacteriia",
                "Saccharimonadia") ~ "Minisyncoccota",
    Class %in% c("Mollicutes", "Izemoplasmatia") ~ "Mycoplasmatota",
    TRUE ~ "Otros"
  ))

library(ComplexHeatmap)
library(circlize)
library(dplyr)

df_plot <- dfc %>% 
  filter(perc_num > 0) %>%
  arrange(Phylum, desc(perc_num))

mat <- as.matrix(df_plot$perc_num)
rownames(mat) <- df_plot$Class  

col_fun = colorRamp2(c(0, max(df_plot$perc_num, na.rm = TRUE)), c("white", "darkgreen"))

row_ha = rowAnnotation(
  "Total Genomes" = anno_barplot(
    df_plot$Total_Genomes, 
    fill = "steelblue", 
    width = unit(3, "cm")
  ),
  show_annotation_name = TRUE
)

ht <- Heatmap(mat, 
              name = "% GAD\nComplete",
              col = col_fun,
              
              row_split = df_plot$Phylum,       
              cluster_rows = FALSE,         
              row_title_rot = 0,            
              row_title_gp = gpar(fontsize = 10, fontface = "bold"),
              row_names_side = "left",
              row_names_gp = gpar(fontface = "italic"), 
              
              right_annotation = row_ha,
              rect_gp = gpar(col = "white", lwd = 1), 
              cell_fun = function(j, i, x, y, w, h, col) {
                grid.text(sprintf("%.1f%%", mat[i, j]), x, y, 
                          gp = gpar(fontsize = 8, 
                                    col = ifelse(mat[i, j] > 20, "white", "black")))
              })

draw(ht)

png("results/reference_plotting/Heatmap_GAD_Taxonomia_Bacteria.png", width = 2800, height = 3500, res = 300)
draw(ht)
dev.off()


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
    log2OR = log2(OR)
  ) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p_val, method = "fdr"))

enrich_plot_data <- enrich_df %>%
  filter(perc_num > 0 & OR > 0) %>% 
  arrange(Phylum, log2OR) %>%
  mutate(Class = factor(Class, levels = Class)) 

enrich_plot <- ggplot(enrich_plot_data, aes(x = Class, y = log2OR, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  geom_text(aes(label = ifelse(fdr < 0.05, "*", "")), 
            hjust = ifelse(enrich_plot_data$log2OR > 0, -0.5, 1.5), 
            vjust = 0.8, size = 6) +
  coord_flip() +
  scale_fill_brewer(palette = "Set3") + 
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Enrichment of GAD System by Taxonomic Class",
    subtitle = "Only classes with presence (>0%) | Bars: Log2 Odds Ratio | Asterisk (*): FDR < 0.05",
    x = "Taxonomic Class",
    y = "Log2 Odds Ratio (Enrichment > 0)",
    fill = "Class"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 9, face = "italic"), 
    legend.position = "right"
  )

print(enrich_plot)
png("results/reference_plotting/Enrichment_GAD_Fisher_Bacteria.png", width = 2400, height = 3000, res = 300)
print(enrich_plot)
dev.off()

df_arch <- read.table("RESUMEN_ARCHAEA_PHYLUMS.tsv", sep = "\t", header=TRUE)


df_arch2 <- df_arch %>% 
  select(Phylum, COMPLETE, PARTIAL_ENZYME, PARTIAL_TRANSPORTER, ABSENT) %>% 
  pivot_longer(cols= -Phylum, names_to ="Status", values_to = "Count")

df_arch2$Status <- factor(df_arch2$Status, 
                         levels = c("COMPLETE", "PARTIAL_ENZYME", "PARTIAL_TRANSPORTER", "ABSENT"))

proporciones2 <- ggplot(df_arch2, aes(x = reorder(Phylum, Count, sum), y = Count, fill = Status)) +
  geom_bar(stat = "identity", position = "fill") + 
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) + 
  scale_fill_manual(values = c("lightgreen", "lightblue", "lightyellow", "salmon")) +
  coord_flip() + 
  labs(title = "GAD system distribution in Archaea",
       subtitle = "Relative proportion",
       x = "Phylum", y = "Proportion (%)", fill = "GAD state") +
  theme_minimal()

absolutos2 <- ggplot(df_arch, aes(x = reorder(Phylum, Total_Genomes), y = Total_Genomes)) +
  geom_bar(stat = "identity", fill = "grey80") +
  geom_text(aes(label = Total_Genomes), hjust = -0.1, size = 3) +
  coord_flip() +
  labs(title = "Genomic representation",
       subtitle = paste("Total:", sum(df_arch$Total_Genomes)),
       x = "", y = "Total genomes") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )

panel_final2 <- proporciones2 + absolutos2 + 
  plot_layout(widths = c(2, 1), guides = "collect") & 
  theme(legend.position = "bottom")

ggsave("RefSeq_plots_archaea.png", panel_final2, width = 20, height = 10)


df_refseq_raw <- read.table("RESUMEN_CLASSES.tsv", header=TRUE, sep="\t")

df_refseq_prep <- df_refseq_raw %>%
  mutate(
    Phylum_Unified = case_when(
      Class %in% c("Bacteroidia", "Cytophagia", "Chitinophagia","Flavobacteriia", 
                   "Saprospiria", "Sphingobacteriia", "Ignavibacteria") ~ "Bacteroidota",
      Class %in% c("Bacilli", "Clostridia", "Culicoidibacteria", "Negativicutes",
                   "Erysipelotrichia", "Limnochordia", "Syntrophomonadia",
                   "Thermaerobacteria", "Thermolithobacteria") ~ "Bacillota / Firmicutes",
      Class %in% c("Alphaproteobacteria", "Betaproteobacteria", 
                   "Gammaproteobacteria", "Deltaproteobacteria", "Epsilonproteobacteria",
                   "Zetaproteobacteria", "Hydrogenophilalia", "Oligoflexia", "Acidithiobacillales", "Magnetococcia") ~ "Pseudomonadota / Proteobacteria",
      Class %in% c("Actinomycetes", "Coriobacteriia", "Acidimicrobiia",
                   "Nitriliruptoria", "Rubrobacteria", "Thermoleophilia") ~ "Actinomycetota",
      Class %in% c("Fusobacteriia") ~ "Fusobacteriota",
      Class %in% c("Chlamydiia") ~ "Chlamydiota",
      Class %in% c("Spirochaetia") ~ "Spirochaetota",
      Class %in% c("Cyanophyceae", "Vampirovibrionophyceae") ~ "Cyanobacteriota",
      Class %in% c("Planctomycetia", "Phycisphaerae") ~ "Planctomycetota",
      Class %in% c("Anaerolineae", "Ardenticatenia", "Caldilineae", "Chloroflexia", "Dehalococcoidia",
                   "Ktedonobacteria", "Tepidiformia", "Thermoflexia") ~ "Chloroflexota",
      Class %in% c("Blastocatellia", "Holophagae","Thermoanaerobaculia", "Vicinamibacteria") ~ "Acidobacteriota",
      Class %in% c("Aquificia", "Desulfurobacteriia", "Thermosulfidibacteria", "Atribacteria") ~ "Aquificota",
      Class %in% c("Balneolia") ~ "Balneolota",
      Class %in% c("Caldisericia") ~ "Caldisericota",
      Class %in% c("Calditrichia") ~ "Calditrichota",
      Class %in% c("Chlorobia") ~ "Chlorobiota",
      Class %in% c("Chrysiogenia") ~ "Chrysiogenota",
      Class %in% c("Coprothermobacteria") ~ "Coprothermobacterota",
      Class %in% c("Deferribacteres") ~ "Deferribacterota",
      Class %in% c("Dictyoglomeria") ~ "Dyctyoglomerota",
      Class %in% c("Elusimicrobia", "Endomicrobiia") ~ "Elusimicrobiota",
      Class %in% c("Chitinispirillia", "Chitinivibrionia", "Fibrobacteria") ~ "Fibrobacterota",
      Class %in% c("Gemmatimonadetes", "Longimicrobia") ~ "Gemmatimonadota",
      Class %in% c("Kiritimatiellae", "Tichowtungiia") ~ "Kiritimatiellota",
      Class %in% c("Oligosphaeria", "Lentisphaeria") ~ "Lentisphaerota",
      Class %in% c("Bradymonadia", "Myxococcia", "Polyangiia") ~ "Myxococcota",
      Class %in% c("Nitrospinia") ~ "Nitrospinota",
      Class %in% c("Nitrospiria", "Thermodesulfovibrionia") ~ "Nitrospirota",
      Class %in% c("Rhodothermia") ~ "Rhodothermota",
      Class %in% c("Deferrisomatia", "Desulfarculia", "Desulfobaccia", "Desulfobacteria",
                   "Desulfobulbia", "Desulfomonilia", "Desulfovibrionia", "Desulfuromonadia",
                   "Dissulfuribacteria", "Syntrophia", "Syntrophobacteria", 
                   "Syntrophorhabdia", "Thermodesulfobacteria") ~ "Thermodesulfobacteriota",
      Class %in% c("Thermomicrobia") ~ "Thermomicrobiota",
      Class %in% c("Opitutia", "Verrucomicrobiia") ~ "Verrucomicrobiota",
      Class %in% c("Deinococci") ~ "Deinococcota",
      Class %in% c("Synergistia") ~ "Synergistota",
      Class %in% c("Thermotogae") ~ "Thermotogota",
      Class %in% c("Armatimonadia", "Chthonomonadia", "Fimbriimonadia") ~ "Armatimonadota",
      Class %in% c("Gracilibacteria", "Microgenomatia", "Minisyncoccia", "Patescibacteriia",
                   "Saccharimonadia") ~ "Minisyncoccota",
      Class %in% c("Mollicutes", "Izemoplasmatia") ~ "Mycoplasmatota",
      TRUE ~ "Otros"
    ),
    Biome = "Reference (RefSeq)"
  ) %>%
  select(Biome, Phylum_Unified, COMPLETE, PARTIAL_ENZYME, PARTIAL_TRANSPORTER, ABSENT) %>%
  pivot_longer(cols = c(COMPLETE, PARTIAL_ENZYME, PARTIAL_TRANSPORTER, ABSENT), 
               names_to = "Status", values_to = "Count") %>%
  uncount(Count) %>%
  select(Biome, Phylum = Phylum_Unified, Status)


df_cfmd_prep_tax <- df_tax_clean %>%
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
      TRUE                                          ~ "Food: Other"
    ),
    Status = System_Status,
    Phylum_Unified = case_when(
      str_detect(Phylum, "Firmicutes|Bacillota") ~ "Bacillota / Firmicutes",
      str_detect(Phylum, "Proteobacteria|Pseudomonadota") ~ "Pseudomonadota / Proteobacteria",
      str_detect(Phylum, "Actinobacteriota|Actinomycetota") ~ "Actinomycetota",
      str_detect(Phylum, "Bacteroidetes|Bacteroidota") ~ "Bacteroidota",
      TRUE ~ as.character(Phylum)
    )
  ) %>%
  select(Biome, Phylum = Phylum_Unified, Status)

df_mgnify_prep_tax <- df_mgnify %>%
  filter(Domain == "Bacteria") %>% 
  mutate(
    Status = case_when(
      Status == "PARTIAL_ENZYME_ONLY" ~ "PARTIAL_ENZYME",
      Status == "PARTIAL_TRANSPORTER_ONLY" ~ "PARTIAL_TRANSPORTER",
      TRUE ~ Status
    ),
    Biome = paste0("MGnify: ", str_to_title(gsub("-", " ", Biome))),
    Phylum_Unified = case_when(
      str_detect(Phylum, "Firmicutes|Bacillota") ~ "Bacillota / Firmicutes",
      str_detect(Phylum, "Proteobacteria|Pseudomonadota") ~ "Pseudomonadota / Proteobacteria",
      str_detect(Phylum, "Actinobacteriota|Actinomycetota") ~ "Actinomycetota",
      str_detect(Phylum, "Bacteroidetes|Bacteroidota") ~ "Bacteroidota",
      TRUE ~ as.character(Phylum)
    )
  ) %>%
  select(Biome, Phylum = Phylum_Unified, Status)


df_triple_unificado <- bind_rows(df_refseq_prep, df_cfmd_prep_tax, df_mgnify_prep_tax) %>%
  filter(!is.na(Biome) & Biome != "") %>%
  filter(!is.na(Phylum) & Phylum != "" & Phylum != "Otros")

df_enrich_counts <- df_triple_unificado %>%
  group_by(Biome, Phylum) %>%
  summarise(
    COMPLETE = sum(Status == "COMPLETE"),
    Total = n(),
    .groups = 'drop'
  ) %>%
  filter(Total >= 5)

t_comp_global <- sum(df_enrich_counts$COMPLETE)
t_gen_global  <- sum(df_enrich_counts$Total)

enrich_triple_data <- df_enrich_counts %>%
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
  mutate(Label_Eje = paste0(Biome, " | ", Phylum)) %>%
  filter(fdr < 0.05 & is.finite(log2OR) & COMPLETE > 0) %>%
  arrange(log2OR) %>%
  mutate(Label_Eje = factor(Label_Eje, levels = Label_Eje))


plot_triple_final <- ggplot(enrich_triple_data, aes(x = Label_Eje, y = log2OR, fill = Phylum)) +
  geom_bar(stat = "identity", color = "black", size = 0.1) +
  geom_text(aes(label = "*"), 
            hjust = ifelse(enrich_triple_data$log2OR > 0, -0.4, 1.4), 
            size = 6, fontface = "bold", vjust = 0.7) +
  coord_flip() +
  scale_fill_brewer(palette = "Set3") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
  labs(
    title = "Macro-Enrichment of GAD System (Joint Bacteria Analysis)", 
    subtitle = "Fisher's Exact Test Combining RefSeq, cFMD, and MGnify Data | Asterisk (*): FDR < 0.05",
    x = "Dataset Environment | Bacterial Phylum", 
    y = "Log2 Odds Ratio (Enrichment > 0)",
    fill = "Phylum"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8, face = "bold"), 
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(plot_triple_final)
ggsave("results/Enrichment_TRIPLE_Global_Bacterias.png", plot_triple_final, width = 15, height = 14, bg = "white")
