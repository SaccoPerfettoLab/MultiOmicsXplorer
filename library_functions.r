library(shiny)
library(shinydashboard)
library(shinyBS)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(fs)
library(ggsignif)
library(shinycssloaders)
library(png)
library(plotly)
library(gridExtra)
library(ggpubr)
library(shinyjs)
library(DT)
library(openxlsx)
library(plumber)
library(httr)
library(readr) 



# color mapping for the plots 
stage_colors <- c("1" = "#66c2a5", "2" = "#fc8d62", "3" = "#8da0cb", "4" = "#e78ac3")

tumor_colors <- c("Brca" = "#e06666", "Ccrcc" = "#f9cb9c", "Coad" = "#ffe599",
                   "Gbm" = "#76a5af", "Hnscc" = "#9fc5e8",
                  "Lscc" = "#8e7cc3", "Luad" = "#d5a6bd", "Ov" = "#d9ead3",
                  "Pdac" = "#fff2cc", "Ucec" = "#b6d7a8") 

normal_color <- "#C0C0C0"  


#### ANALYTE ABUNDANCE ####
load_data <- function(data_type, gene) {
  data_type_uc <- tools::toTitleCase(data_type)
  data_type_prefix <- switch(data_type_uc,
                             "Proteomics" = "Prot",
                             "Phosphoproteomics" = "Phospho",
                             "Transcriptomics" = "Transc")
  
  genes_file <- paste0(data_type_prefix, "_genes.tsv")
  if (!fs::file_exists(genes_file)) {
    cat("Error: Gene list file not found.\n")
    return(NULL)
  }
  
  gene_list <- readr::read_tsv(genes_file)$gene_name
  
  gene_letter <- substr(gene, 1, 1) ## it takes the first letter of the name of the analyte cause files 
                                    ## are divided per first letter
  
  file_name <- paste0("Gene_expression/", data_type_prefix, "_data_", gene_letter, ".tsv") ## Gene_expression is the folder containing all 
                                                                                           ## the files of analyte abundance 
  
  cat("File path:", file_name, "\n")
  if (fs::file_exists(file_name)) {
    data <- readr::read_tsv(file_name)
    return(data)
  } else {
    return(NULL)
  }
}

#### ABUNDANCE PLOT GENERATION ####

# all the combination are managed to have: tum vs tum with >=1 stages; tum in >1 stages; tum vs norm in all (average) or >1 stages
generate_plot_abundance <- function(data, gene, tumors, stages, data_type, 
                                    comparison_type, analyte_option = "all", 
                                    show_points = FALSE, show_significance = FALSE) {
  

  if (comparison_type == "Tumor vs Tumor") {
    
    if ("all" %in% stages) {
      filtered_data <- data %>%
        filter(Name == gene, Tumor_Type %in% tumors, !grepl("\\.N$", Patient_ID)) %>%
        group_by(Tumor_Type, Patient_ID) %>%
        summarise(z_score = mean(z_score, na.rm = TRUE), .groups = 'drop')
      
      if (analyte_option == "significant") {
        filtered_data <- filtered_data %>% filter(z_score < -1.96 | z_score > 1.96)
        if (nrow(filtered_data) == 0) {
          message("No significant analytes available for Tumor vs Tumor comparison.")
          return(NULL)
        }
      }
      
      available_stages <- data %>%
        filter(Name == gene, Tumor_Type %in% tumors, !grepl("\\.N$", Patient_ID)) %>%
        pull(Stage) %>%
        unique() %>%
        sort()
      
      title <- paste(
        "Analyte:", gene,
        "-", "Tumors:", paste(tumors, collapse = ", "),
        "- Stages: All (", paste(available_stages, collapse = ", "), ")"
      )
      
      p <- ggplot(filtered_data, aes(x = Tumor_Type, y = z_score, fill = Tumor_Type)) +
        geom_boxplot() +
        labs(title = title, x = "Tumor Type", y = "Abundance") +
        scale_fill_manual(values = tumor_colors) +
        theme_minimal() +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.line = element_line(color = "gray"),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 12, face = "bold")
        )
      
      if (show_points) {
        p <- p + 
          geom_point(aes(fill = Tumor_Type), 
                     position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.8, 
                     shape = 21, color = "#000000", stroke = 0.5)
      }
      
      #### signif tum vs tum ####
      if (show_significance) {
        valid_tumors <- unique(filtered_data$Tumor_Type[filtered_data$Tumor_Type %in% tumors])
        
        # Test for significance by tumor type
        if (length(valid_tumors) >= 2) {
          comparisons <- combn(valid_tumors, 2, simplify = FALSE)
          p <- p + stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif")
        }
      } 
    } else {
      filtered_data <- data %>%
        filter(Name == gene, Tumor_Type %in% tumors, Stage %in% stages, !grepl("\\.N$", Patient_ID))
      
      if (nrow(filtered_data) == 0) {
        message("No data available for Tumor vs Tumor comparison.")
        return(NULL)
      }
      
      if (analyte_option == "significant") {
        filtered_data <- filtered_data %>% filter(z_score < -1.96 | z_score > 1.96)
        if (nrow(filtered_data) == 0) {
          message("No significant analytes available for Tumor vs Tumor comparison.")
          return(NULL)
        }
      }
      
      comparisons_by_stage <- filtered_data %>%
        group_by(Stage) %>%
        summarise(comparisons = list(if (length(unique(Tumor_Type)) >= 2) combn(unique(Tumor_Type), 2, simplify = FALSE) else list()), .groups = 'drop')
      
      title <- paste("Analyte:", gene, "-", "Tumors:", paste(tumors, collapse = ", "), "- Stages:", paste(stages, collapse = ", "))
      p <- ggplot(filtered_data, aes(x = Tumor_Type, y = z_score, fill = Tumor_Type)) +
        geom_boxplot() +
        labs(title = title, x = "Tumor Type", y = "Abundance") +
        scale_fill_manual(values = tumor_colors) +
        theme_minimal() +
        facet_wrap(~ Stage, nrow = 1) +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.line = element_line(color = "gray"),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 12, face = "bold")
        )
      
      if (show_points) {
        p <- p + 
          geom_point(aes(fill = Tumor_Type), 
                     position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.8, 
                     shape = 21, color = "#000000", stroke = 0.5)
      }
      
      if (show_significance) {
        valid_tumors <- unique(filtered_data$Tumor_Type[filtered_data$Tumor_Type %in% tumors])
        unique_stages <- unique(filtered_data$Stage)
        
        # One tumor - more stages
        if (length(valid_tumors) == 1 && length(unique_stages) > 1) {
          filtered_data$Tumor_Stage <- paste(filtered_data$Tumor_Type, filtered_data$Stage, sep = " - ")
          
          p <- ggplot(filtered_data, aes(x = Tumor_Stage, y = z_score, fill = Tumor_Type)) +
            geom_boxplot() +
            labs(title = title, x = "Tumor Type - Stage", y = "Abundance") +
            scale_fill_manual(values = tumor_colors) +
            theme_minimal() +
            theme(
              legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
              axis.text.y = element_text(size = 14),
              axis.line = element_line(color = "gray"),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 12, face = "bold")
            )
          
          comparisons <- combn(unique_stages, 2, simplify = FALSE)
          
          valid_comparisons <- comparisons[sapply(comparisons, function(x) {
            sum(filtered_data$Stage == x[1]) > 0 & sum(filtered_data$Stage == x[2]) > 0
          })]
          
          if (length(valid_comparisons) > 0) {
            p <- p + stat_compare_means(method = "wilcox.test", 
                                        comparisons = valid_comparisons, 
                                        aes(group = Stage),  
                                        label = "p.signif",  
                                        size = 5)
          }
        } 
        else if (length(valid_tumors) >= 2) {
          # more tumors - comparison for each stage
          p <- ggplot(filtered_data, aes(x = Tumor_Type, y = z_score, fill = Tumor_Type)) +
            geom_boxplot() +
            facet_wrap(~Stage, scales = "fixed", nrow = 1) +
            labs(title = title, x = "Tumor Type", y = "Abundance") +
            scale_fill_manual(values = tumor_colors) +
            theme_minimal() +
            theme(
              legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
              axis.text.y = element_text(size = 14),
              axis.line = element_line(color = "gray"),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 12, face = "bold")
            )
          
          for (stage in unique_stages) {
            stage_data <- filtered_data[filtered_data$Stage == stage, ]
            stage_tumors <- unique(stage_data$Tumor_Type)
            
            if (length(stage_tumors) >= 2) {
              comparisons <- combn(stage_tumors, 2, simplify = FALSE)
              
              valid_comparisons <- comparisons[sapply(comparisons, function(x) {
                sum(stage_data$Tumor_Type == x[1]) > 0 & sum(stage_data$Tumor_Type == x[2]) > 0
              })]
              
              if (length(valid_comparisons) > 0) {
                p <- p + stat_compare_means(method = "wilcox.test", 
                                            comparisons = valid_comparisons, 
                                            aes(group = Tumor_Type),  
                                            label = "p.signif",  
                                            size = 5, 
                                            data = stage_data)
              }
            }
          }
        }
      }
      
        if (show_points) {
          p <- p + 
            geom_point(aes(fill = Tumor_Type), 
                       position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.8, 
                       shape = 21, color = "#000000", stroke = 0.5)
        }
        
      }
    
  } else if (comparison_type == "Tumor vs Normal") {
    
    if ("all" %in% stages) {
      filtered_data <- data %>%
        filter(Name == gene, Tumor_Type %in% tumors) %>%
        group_by(Tumor_Type, Patient_ID) %>%
        summarise(z_score = mean(z_score, na.rm = TRUE), .groups = 'drop')
      
      if (nrow(filtered_data) == 0) {
        message("No data available for Tumor vs Normal comparison.")
        return(NULL)
      }
      
      if (analyte_option == "significant") {
        filtered_data <- filtered_data %>% filter(z_score < -1.96 | z_score > 1.96)
        if (nrow(filtered_data) == 0) {
          message("No significant analytes available for Tumor vs Normal comparison.")
          return(NULL)
        }
      }
      
      filtered_data <- filtered_data %>%
        mutate(Sample_Type = ifelse(grepl("\\.N$", Patient_ID), "Normal", Tumor_Type))
      
      available_stages <- data %>%
        filter(Name == gene, Tumor_Type %in% tumors, !grepl("\\.N$", Patient_ID)) %>%
        pull(Stage) %>%
        unique() %>%
        na.omit() %>%
        sort()
      
      title <- paste0("Analyte: ", gene, 
                      " – Tumor: ", paste(tumors, collapse = ", "), 
                      " – Stages: All (", paste(available_stages, collapse = ", "), ")")
      
      p <- ggplot(filtered_data, aes(x = Sample_Type, y = z_score, fill = Sample_Type)) +
        geom_boxplot() +
        labs(title = title, x = "Sample Type", y = "Abundance") +
        scale_fill_manual(values = c("Normal" = "gray", tumor_colors)) +
        theme_minimal() +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.line = element_line(color = "gray"),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 12, face = "bold")
        )
      
      if (show_points) {
        p <- p + 
          geom_point(aes(fill = Sample_Type), 
                     position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.8, 
                     shape = 21, color = "#000000", stroke = 0.5)
      }
      
      if (show_significance) {
        valid_tumors <- unique(filtered_data$Tumor_Type)
        
        comparisons <- lapply(valid_tumors, function(tumor) c("Normal", tumor))
        
        if (length(comparisons) > 0) {
          p <- p + stat_compare_means(
            method = "wilcox.test",
            comparisons = comparisons,
            label = "p.signif"
          )
        }
      }
    } else {
      filtered_data <- data %>%
        filter(Name == gene, Tumor_Type %in% tumors, Stage %in% stages | grepl("\\.N$", Patient_ID))
      
      if (nrow(filtered_data) == 0) {
        message("No data available for Tumor vs Normal comparison.")
        return(NULL)
      }
      
      if (analyte_option == "significant") {
        filtered_data <- filtered_data %>% filter(z_score < -1.96 | z_score > 1.96)
        if (nrow(filtered_data) == 0) {
          message("No significant analytes available for Tumor vs Normal comparison.")
          return(NULL)
        }
      }
      
      filtered_data <- filtered_data %>%
        mutate(Sample_Type = ifelse(grepl("\\.N$", Patient_ID), "Normal", paste(Tumor_Type, Stage, sep = " - ")))
      
      title <- paste("Analyte:", gene, "- Tumor:", paste(tumors, collapse = ", "), "- Stages:", paste(stages, collapse = ", "))
      p <- ggplot(filtered_data, aes(x = Sample_Type, y = z_score, fill = ifelse(Sample_Type == "Normal", "Normal", Tumor_Type))) +
        geom_boxplot() +
        labs(title = title, x = "Sample Type", y = "Abundance") +
        scale_fill_manual(values = c("Normal" = "gray", tumor_colors)) +
        theme_minimal() +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.line = element_line(color = "gray"),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 12, face = "bold")
        )
      
      if (show_points) {
        p <- p + 
          geom_point(aes(fill = Sample_Type), 
                     position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.8, 
                     shape = 21, color = "#000000", stroke = 0.5)
      }
      
      if (show_significance) {
        valid_sample_types <- unique(filtered_data$Sample_Type)
        
        if (length(valid_sample_types) < 2) {
          message("Not enough valid sample types for significance testing.")
          return(NULL)
        }
        
        comparisons <- lapply(valid_sample_types[valid_sample_types != "Normal"], function(sample) c("Normal", sample))
        
        if (length(comparisons) > 0) {
          p <- p + stat_compare_means(
            method = "wilcox.test",
            comparisons = comparisons,
            label = "p.signif"
          )
        } else {
          message("There are no valid tumor comparisons for significance testing.")
        }
      }
      
    }
  }
  return(p)
}


#### PROTEIN ACTIVITY ####
load_protein_choices <- function() {
  gene_names_file <- "Prot_Act_names.tsv"
  
  if (!file.exists(gene_names_file)) {
    cat("Error: Gene names file not found.\n")
    return(NULL)
  }
  
  gene_names <- readr::read_tsv(gene_names_file)$gene_name
  return(gene_names)
}

#### load data protein ####
load_data_protein <- function(protein) {
  gene_names_file <- "Prot_Act_names.tsv"
  
  if (!file.exists(gene_names_file)) {
    cat("Error: Protein names file not found.\n")
    return(NULL)
  }
  
  gene_names <- readr::read_tsv(gene_names_file)$gene_name
  
  if (!protein %in% gene_names) {
    cat("Error: Protein not found in the gene names list.\n")
    return(NULL)
  }
  
  gene_letter <- substr(protein, 1, 1)
  
  file_name <- paste0("Protein_activity/Protein_activity_", gene_letter, ".tsv")
  
  cat("File path:", file_name, "\n")
  if (file.exists(file_name)) {
    data <- readr::read_tsv(file_name)
    return(data)
  } else {
    cat("Error: Protein Activity data file not found.\n")
    return(NULL)
  }
}

#### PROTEIN ACTIVITY GENERATION PLOT ####

# all the combination are managed to have: tum vs tum with >=1 stages; tum in >1 stages; tum vs norm in all (average) or >1 stages
generate_plot_protein_activity <- function(data, protein, tumors, stages, 
                                           comparison_type_prot, 
                                           show_points_protein = FALSE, 
                                           show_significance = FALSE) {
  
  if (comparison_type_prot == "Tumor vs Tumor") {.  
    
    if ("all" %in% stages) {
      filtered_data <- data %>%
        filter(Name == protein, Tumor_Type %in% tumors, !grepl("\\.N$", Patient_ID)) %>%
        group_by(Tumor_Type, Patient_ID) %>%
        summarise(predicted_activity = mean(predicted_activity, na.rm = TRUE), .groups = 'drop')
      
      head(filtered_data)
      if (nrow(filtered_data) == 0) {
        message("No data available for the selected protein and tumors with 'all' stages.")
        return(NULL)
      }
      
      available_stages <- data %>%
        filter(Name == protein, Tumor_Type %in% tumors, !grepl("\\.N$", Patient_ID)) %>%
        pull(Stage) %>%
        unique() %>%
        na.omit() %>%
        sort()
      
      title <- paste0("Protein: ", protein, 
                      " – Tumor: ", paste(tumors, collapse = ", "), 
                      " – Stages: All (", paste(available_stages, collapse = ", "), ")")    
      
      p <- ggplot(filtered_data, aes(x = Tumor_Type, y = predicted_activity, fill = Tumor_Type)) +
        geom_boxplot() +
        labs(title = title, x = "Tumor Type", y = "Predicted activity") +
        scale_fill_manual(values = tumor_colors) +
        theme_minimal() +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.line = element_line(color = "gray"),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 12, face = "bold")
        )
      
      ## here we have the wilcox test between tumors in the all option
      if (show_significance && length(tumors) >= 2) {
        valid_tumors <- tumors[sapply(tumors, function(tumor) any(filtered_data$Tumor_Type == tumor))]
        
        if (length(valid_tumors) >= 2) {
          p <- p + stat_compare_means(method = "wilcox.test", 
                                      comparisons = combn(valid_tumors, 2, simplify = FALSE), 
                                      label = "p.signif")
        } else {
          message("There are not enough valid tumors for significance testing.")
        }
      }
      
      if (show_points_protein) {
        p <- p + 
          geom_point(aes(fill = Tumor_Type), 
                     position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.8, 
                     shape = 21, color = "#000000", stroke = 0.5)
      } 
      
    } else {
      filtered_data <- data %>%
        filter(Name == protein, Tumor_Type %in% tumors, Stage %in% stages, !grepl("\\.N$", Patient_ID))
      
      if (nrow(filtered_data) == 0) {
        message("No data available for the selected protein, tumors, and stages.")
        return(NULL)
      }
      
      title <- paste("Protein:", protein, "-", "Tumors:", paste(tumors, collapse = ", "), "- Stages:", paste(stages, collapse = ", "))
      
      p <- ggplot(filtered_data, aes(x = Tumor_Type, y = predicted_activity, fill = Tumor_Type)) +
        geom_boxplot() +
        labs(title = title, x = "Tumor Type", y = "Predicted activity") +
        scale_fill_manual(values = tumor_colors) +
        theme_minimal() +
        facet_wrap(~ Stage, nrow = 1) +  
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.line = element_line(color = "gray"),
          axis.title = element_text(size = 12),
          plot.title = element_text(size = 12, face = "bold")
        )
      
      if (show_significance) {
        valid_tumors <- unique(filtered_data$Tumor_Type[filtered_data$Tumor_Type %in% tumors])
        unique_stages <- unique(filtered_data$Stage)
        
        ## comparison between stages if one tumor selected
        if (length(valid_tumors) == 1 && length(unique_stages) > 1) { # es brca 1-2-3-4
          filtered_data$Tumor_Stage <- paste(filtered_data$Tumor_Type, filtered_data$Stage, sep = " - ")
          
          p <- ggplot(filtered_data, aes(x = Tumor_Stage, y = predicted_activity, fill = Tumor_Type)) +
            geom_boxplot() +
            labs(title = title, x = "Tumor Type - Stage", y = "Predicted activity") +
            scale_fill_manual(values = tumor_colors) +
            theme_minimal() +
            theme(
              legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
              axis.text.y = element_text(size = 14),
              axis.line = element_line(color = "gray"),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 12, face = "bold")
            )
          
          comparisons <- combn(unique(filtered_data$Tumor_Stage), 2, simplify = FALSE)
          
          valid_comparisons <- comparisons[sapply(comparisons, function(x) {
            sum(filtered_data$Tumor_Stage == x[1]) > 0 | sum(filtered_data$Tumor_Stage == x[2]) > 0
          })]
          
          ## wilcox for the same tumor in more stages, es brca-1 vs brca-2 / brca-1 vs brca-3 etc
          if (length(valid_comparisons) > 0) {
            p <- p + stat_compare_means(method = "wilcox.test", 
                                        comparisons = valid_comparisons, 
                                        aes(group = Tumor_Stage),  
                                        label = "p.signif",  
                                        size = 5)
          }
        } 
        else if (length(valid_tumors) >= 2) {
          
          ## comparison with tumors and stages, es Brca, Ccrcc, Coad and stages 1,2,3
          p <- ggplot(filtered_data, aes(x = Tumor_Type, y = predicted_activity, fill = Tumor_Type)) +
            geom_boxplot() +
            facet_wrap(~Stage, scales = "fixed", nrow = 1) +
            labs(title = title, x = "Tumor Type", y = "Predicted activity") +
            scale_fill_manual(values = tumor_colors) +
            theme_minimal() +
            theme(
              legend.position = "none",
              axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
              axis.text.y = element_text(size = 14),
              axis.line = element_line(color = "gray"),
              axis.title = element_text(size = 14),
              plot.title = element_text(size = 12, face = "bold")
            )
          
          for (stage in unique_stages) {
            stage_data <- filtered_data[filtered_data$Stage == stage, ]
            stage_tumors <- unique(stage_data$Tumor_Type)
            
            if (length(stage_tumors) >= 2) {
              comparisons <- combn(stage_tumors, 2, simplify = FALSE)
              
              valid_comparisons <- comparisons[sapply(comparisons, function(x) {
                sum(stage_data$Tumor_Type == x[1]) > 0 & sum(stage_data$Tumor_Type == x[2]) > 0
              })]
               
            ## wilcox between tumors in each stage selected
              if (length(valid_comparisons) > 0) {
                p <- p + stat_compare_means(method = "wilcox.test", 
                                            comparisons = valid_comparisons, 
                                            aes(group = Tumor_Type),  
                                            label = "p.signif",  
                                            size = 5, 
                                            data = stage_data)
              }
            }
          }
        }
      }
      
      if (show_points_protein) {
        p <- p + 
          geom_point(aes(fill = Tumor_Type), 
                     position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.8, 
                     shape = 21, color = "#000000", stroke = 0.5)
      } 
    } # end "Tumor vs Tumor" block -> start with "Tumor vs Normal"
    
  }else if (comparison_type_prot== "Tumor vs Normal") {
    
    if ("all" %in% stages) {
      filtered_data <- data %>%
        filter(Name == protein, Tumor_Type %in% tumors) %>%
        group_by(Tumor_Type, Patient_ID) %>%
        summarise(predicted_activity = mean(predicted_activity, na.rm = TRUE), .groups = 'drop')
      
      if (nrow(filtered_data) == 0) {
        message("No data available for Tumor vs Normal comparison.")
        return(NULL)
      }
      
      filtered_data <- filtered_data %>%
        mutate(Sample_Type = ifelse(grepl("\\.N$", Patient_ID), "Normal", Tumor_Type))
      
      available_stages <- data %>%
        filter(Name == protein, Tumor_Type %in% tumors, !grepl("\\.N$", Patient_ID)) %>%
        pull(Stage) %>%
        unique() %>%
        na.omit() %>%
        sort()
      
      title <- paste0("Protein: ", protein, 
                      " – Tumor: ", paste(tumors, collapse = ", "), 
                      " – Stages: All (", paste(available_stages, collapse = ", "), ")")    
      ## all stages for the 1 tumor vs normal
      p <- ggplot(filtered_data, aes(x = Sample_Type, y = predicted_activity, fill = Sample_Type)) +
        geom_boxplot() +
        labs(title = title, x = "Sample Type", y = "Predicted activity") +
        scale_fill_manual(values = c("Normal" = "gray", tumor_colors)) +
        theme_minimal() +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.line = element_line(color = "gray"),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 12, face = "bold")
        )
      
      if (show_points_protein) {
        p <- p + 
          geom_point(aes(fill = Sample_Type), 
                     position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.8, 
                     shape = 21, color = "#000000", stroke = 0.5)
      }
      
      if (show_significance) {
        valid_tumors <- unique(filtered_data$Tumor_Type)
         ## wilcox all vs normal
        comparisons <- lapply(valid_tumors, function(tumor) c("Normal", tumor))
        
        if (length(comparisons) > 0) {
          p <- p + stat_compare_means(
            method = "wilcox.test",
            comparisons = comparisons,
            label = "p.signif"
          )
        }
      }
    } else {
      filtered_data <- data %>%
        filter(Name == protein, Tumor_Type %in% tumors, Stage %in% stages | grepl("\\.N$", Patient_ID))
      
      if (nrow(filtered_data) == 0) {
        message("No data available for Tumor vs Normal comparison.")
        return(NULL)
      }
      
      filtered_data <- filtered_data %>%
        mutate(Sample_Type = ifelse(grepl("\\.N$", Patient_ID), "Normal", paste(Tumor_Type, Stage, sep = " - ")))
      
      title <- paste("Analyte:", protein, "- Tumor:", paste(tumors, collapse = ", "), "- Stages:", paste(stages, collapse = ", "))
      p <- ggplot(filtered_data, aes(x = Sample_Type, y = predicted_activity, fill = ifelse(Sample_Type == "Normal", "Normal", Tumor_Type))) +
        geom_boxplot() +
        labs(title = title, x = "Sample Type", y = "Predicted activity") +
        scale_fill_manual(values = c("Normal" = "gray", tumor_colors)) +
        theme_minimal() +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.line = element_line(color = "gray"),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 12, face = "bold")
        )
      
      if (show_points_protein) {
        p <- p + 
          geom_point(aes(fill = Sample_Type), 
                     position = position_jitter(width = 0.2, height = 0), size = 1.5, alpha = 0.8, 
                     shape = 21, color = "#000000", stroke = 0.5)
      }
      
      if (show_significance) {
        valid_sample_types <- unique(filtered_data$Sample_Type)
        
        if (length(valid_sample_types) < 2) {
          message("Not enough valid sample types for significance testing.") # if there's not at least 2 valid sample the comparison is not possible
          return(NULL)
        }
        
        comparisons <- lapply(valid_sample_types[valid_sample_types != "Normal"], function(sample) c("Normal", sample))
        
        if (length(comparisons) > 0) {
          p <- p + stat_compare_means(
            method = "wilcox.test",
            comparisons = comparisons,
            label = "p.signif"
          )
        } else {
          message("There are no valid tumor comparisons for significance testing.")
        }
      }
      
    }
  }
  
  return(p)
} 



#### TABLES GENERATION ####

#### Abund norm ####
generate_table_norm_ab <- function(data, gene, tumor_type, stage, analyte_option = "all") {
  tumor_data <- data %>%
    filter(Name == gene) %>%
    mutate(Type = if_else(grepl("\\.N$", Patient_ID), "Normal", "Tumor"))  
  
  tumor_data <- tumor_data %>%
    filter(
      Tumor_Type %in% as.character(tumor_type), 
      Type == "Normal" | (Type == "Tumor" & (is.na(Stage) | "all" %in% stage | Stage %in% stage))
    )
  
  if (analyte_option == "significant") {
    tumor_data <- tumor_data %>% filter(z_score < -1.96 | z_score > 1.96)
  }
  
  return(tumor_data)
}

#### Prot act norm ####
generate_table_norm_ac <- function(data, protein, tumor_type, stage = "all") {
  tumor_data <- data %>%
    filter(Name == protein) %>%
    mutate(Type = if_else(grepl("\\.N$", Patient_ID), "Normal", "Tumor"))

  tumor_data <- tumor_data %>%
    filter(
      Tumor_Type %in% as.character(tumor_type),  
      Type == "Normal" | (Type == "Tumor" & (is.na(Stage) | "all" %in% stage | Stage %in% stage))
    )
  
  return(tumor_data)
}





#### QUERY TO INTERROGATE SIGNOR AND TO CREATE THE LINK FOR PHOSPHOSIGNOR ####
 
query_ptm_residue <- function(uniprot_id, protein_name, residue) {
  base_url <- "https://signor.uniroma2.it/PhosphoSIGNOR/apis/residueSearch.php"
  
  response <- GET(base_url, query = list(id = uniprot_id, residue = residue)) 
  
  if (http_error(response)) {
    cat("Error: Unable to connect to SIGNOR API.\n")
    return(NULL)
  }
  
  result <- content(response, as = "text", encoding = "UTF-8")
  print(result) 
  
  if (grepl("No result found!", result)) {
    cat("Residue not found in SIGNOR database.\n")
    return(NULL)
  }
  
  formatted_id <- paste0(uniprot_id, "_", protein_name, "_", residue)
  
  network_url <- paste0("https://signor.uniroma2.it/PhosphoSIGNOR/results/entity.php?role=residue&ID=", formatted_id) # remember to check sometimes if the name/link of the resource has been changed
  
  return(network_url)
}

