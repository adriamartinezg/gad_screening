#!/usr/bin/env Rscript

library(dplyr)
library(readr)
library(taxize)
library(rentrez)
library(purrr)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("No se encontraron archivos con el patrón GAD_*")

message(paste("Archivos encontrados:", length(args)))

gad_results <- args %>%
  map_df(~read_delim(.x, delim="\t", show_col_types = FALSE)) %>%
  filter(System_Status %in% c("ABSENT", "COMPLETE"))

if (nrow(gad_results) == 0) stop("No hay genomas con el sistema GAD completo.")

unique_accessions <- gsub("\\..*$", "", unique(gad_results$Genome_ID))

gcf_to_taxonomy <- function(gcf) {
  tryCatch({
    res <- entrez_search(db="assembly", term=gcf)
    if(length(res$ids) == 0) return(NULL)
    
    sum <- entrez_summary(db="assembly", id=res$ids[1])
    tx <- classification(sum$taxid, db="ncbi")[[1]]
    
    get_rank <- function(rank_name) {
      val <- tx$name[tx$rank == rank_name]
      return(if(length(val) == 0) NA else val)
    }
    
    return(data.frame(
      Genome_ID_clean = gcf,
      Species = get_rank("species"),
      Phylum  = get_rank("phylum"),
      Class   = get_rank("class"),
      Order   = get_rank("order"),
      Family  = get_rank("family"),
      Genus   = get_rank("genus"),
      stringsAsFactors = FALSE
    ))
  }, error = function(e) return(NULL))
}

message("Consultando taxonomía en NCBI...")
taxonomy_list <- lapply(unique_accessions, gcf_to_taxonomy) %>% bind_rows()

final_df <- gad_results %>%
  mutate(Genome_ID_clean = gsub("\\..*$", "", Genome_ID)) %>%
  left_join(taxonomy_list, by = "Genome_ID_clean") %>%
  select(-Genome_ID_clean)

# Si solo se quieren completos
#  complete_df <- final_df %>% 
#  filter(System_Status == "COMPLETE")

write.table(complete_df, 
            file = "RefSeq_GAD_complete.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)

# Preparación para GTDBTk
#message("Asignando un negativo a cada genoma")

#df_complete_unique <- final_df %>%
#  filter(System_Status == "COMPLETE") %>%
#  group_by(Species) %>%
#  slice(1) %>%  # Nos quedamos con el mejor de cada especie
#  ungroup()

# Identificar qué GÉNEROS tienen al menos un COMPLETE
#generos_con_complete <- unique(df_complete_unique$Genus)

# ABSENT que pertenezcan a los géneros de la lista anterior
#df_absent_controls <- final_df %>%
#  filter(System_Status == "ABSENT") %>%
 # filter(Genus %in% generos_con_complete) %>%
 # group_by(Genus) %>%
 # slice(1) %>%  
 # ungroup()

# Lista para GTDB-Tk
#df_tree_final <- bind_rows(df_complete_unique, df_absent_controls)

#write_delim(df_tree_final, "RS_Final_Tree_Selection_200_1000.tsv", delim="\t")
#write.table(df_tree_final$Genome_ID, "RS_final_ids_for_gtdbtk.txt", 
 #           row.names = FALSE, col.names = FALSE, quote = FALSE)

#message(paste("Selección final lista:", nrow(df_tree_final), "genomas (", 
 #             nrow(df_complete_unique), "Completos y", 
  #            nrow(df_absent_controls), "Controles )"))
