
library(shiny)
library(shinydashboard)
library(shinyBS)
library(shinyjs)
library(shinycssloaders)

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(plotly)
library(gridExtra)

library(DT)
library(openxlsx)
library(plumber)

library(fs)
library(png)

library(httr)
library(readr)

library(SignalingProfiler)

# ----------------------------
# Global color palettes
# ----------------------------

# Stage color mapping (used in some plots / UI legends)
stage_colors <- c("1" = "#66c2a5", "2" = "#fc8d62", "3" = "#8da0cb", "4" = "#e78ac3")

# Tumor color mapping
tumor_colors <- c(
  "Brca"  = "#e06666",
  "Ccrcc" = "#f9cb9c",
  "Coad"  = "#ffe599",
  "Gbm"   = "#76a5af",
  "Hnscc" = "#9fc5e8",
  "Lscc"  = "#8e7cc3",
  "Luad"  = "#d5a6bd",
  "Ov"    = "#d9ead3",
  "Pdac"  = "#fff2cc",
  "Ucec"  = "#b6d7a8"
)

# Default color for non-tumor samples
normal_color <- "#C0C0C0"


# =============================================================
# 1) ANALYTE ABUNDANCE: load data
# =============================================================

# Load analyte abundance table for a given data type and gene.
#
# data_type: one of c("Proteomics", "Phosphoproteomics", "Transcriptomics")
# gene:      analyte name (string)
#
# Returns: a tibble or NULL if not found.
load_data <- function(data_type, gene) {

  # Normalize the user selection (e.g., "proteomics" -> "Proteomics")
  data_type_uc <- tools::toTitleCase(data_type)

  # Prefix used in your repository naming conventions
  data_type_prefix <- switch(
    data_type_uc,
    "Proteomics"        = "Prot",
    "Phosphoproteomics" = "Phospho",
    "Transcriptomics"   = "Transc",
    stop("Unsupported data_type: ", data_type)
  )

  # Gene list file (used mainly for validation / choices)
  genes_file <- paste0(data_type_prefix, "_genes.tsv")
  if (!fs::file_exists(genes_file)) {
    cat("Error: Gene list file not found: ", genes_file, "\n")
    return(NULL)
  }

  # Gene list is read for completeness (even if not used downstream here)
  gene_list <- readr::read_tsv(genes_file, show_col_types = FALSE)$gene_name
  # NOTE: gene_list is currently not used, but you may keep it to validate `gene` if desired.

  # Files are split by the first letter of the analyte name
  gene_letter <- substr(gene, 1, 1)

  # Expected abundance file path
  file_name <- paste0("Gene_expression/", data_type_prefix, "_data_", gene_letter, ".tsv")
  cat("File path:", file_name, "\n")

  if (fs::file_exists(file_name)) {
    return(readr::read_tsv(file_name, show_col_types = FALSE))
  }

  return(NULL)
}


# =============================================================
# 2) ANALYTE ABUNDANCE: plot generation
# =============================================================

# Generate abundance plot (Z-score) for a selected analyte.
#
# comparison_type:
#   - "Tumor vs Tumor"
#   - "Tumor vs Non-tumor"  (Non-tumor == Patient_ID ending with ".N")
#
# analyte_option:
#   - "all"
#   - "significant"  (keeps z_score outside ±1.96)
#
# show_points:       overlay jitter points
# show_significance: add Wilcoxon test p-values (ggpubr::stat_compare_means)
generate_plot_abundance <- function(
    data,
    gene,
    tumors,
    stages,
    data_type,
    comparison_type,
    analyte_option = "all",
    show_points = FALSE,
    show_significance = FALSE
) {

  if (is.null(data) || nrow(data) == 0) {
    message("No data provided to generate_plot_abundance().")
    return(NULL)
  }

  # ----------------------------
  # Tumor vs Tumor
  # ----------------------------
  if (comparison_type == "Tumor vs Tumor") {

    # Case: "all" stages -> collapse per Patient within each tumor
    if ("all" %in% stages) {

      filtered_data <- data %>%
        dplyr::filter(.data$Name == gene, .data$Tumor_Type %in% tumors, !grepl("\\.N$", .data$Patient_ID)) %>%
        dplyr::group_by(.data$Tumor_Type, .data$Patient_ID) %>%
        dplyr::summarise(z_score = mean(.data$z_score, na.rm = TRUE), .groups = "drop")

      if (analyte_option == "significant") {
        filtered_data <- filtered_data %>% dplyr::filter(.data$z_score < -1.96 | .data$z_score > 1.96)
        if (nrow(filtered_data) == 0) {
          message("No significant analytes available for Tumor vs Tumor comparison.")
          return(NULL)
        }
      }

      # Collect stages present in the original (non-normal) subset for title clarity
      available_stages <- data %>%
        dplyr::filter(.data$Name == gene, .data$Tumor_Type %in% tumors, !grepl("\\.N$", .data$Patient_ID)) %>%
        dplyr::pull(.data$Stage) %>%
        unique() %>%
        sort()

      title <- paste(
        "Analyte:", gene,
        "-", "Tumors:", paste(tumors, collapse = ", "),
        "- Stages: All (", paste(available_stages, collapse = ", "), ")"
      )

      p <- ggplot(filtered_data, aes(x = .data$Tumor_Type, y = .data$z_score, fill = .data$Tumor_Type)) +
        geom_boxplot() +
        labs(title = title, x = "Tumor Type", y = "Abundance (Z-Score)") +
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

      if (isTRUE(show_points)) {
        p <- p +
          geom_point(
            aes(fill = .data$Tumor_Type),
            position = position_jitter(width = 0.2, height = 0),
            size = 1.5, alpha = 0.8,
            shape = 21, color = "#000000", stroke = 0.5
          )
      }

      if (isTRUE(show_significance)) {
        valid_tumors <- unique(filtered_data$Tumor_Type[filtered_data$Tumor_Type %in% tumors])
        if (length(valid_tumors) >= 2) {
          comparisons <- combn(valid_tumors, 2, simplify = FALSE)
          p <- p + stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif")
        }
      }

      return(p)
    }

    # Case: specific stages -> keep Stage facet and compare within each stage
    filtered_data <- data %>%
      dplyr::filter(.data$Name == gene, .data$Tumor_Type %in% tumors, .data$Stage %in% stages, !grepl("\\.N$", .data$Patient_ID))

    if (nrow(filtered_data) == 0) {
      message("No data available for Tumor vs Tumor comparison.")
      return(NULL)
    }

    if (analyte_option == "significant") {
      filtered_data <- filtered_data %>% dplyr::filter(.data$z_score < -1.96 | .data$z_score > 1.96)
      if (nrow(filtered_data) == 0) {
        message("No significant analytes available for Tumor vs Tumor comparison.")
        return(NULL)
      }
    }

    title <- paste(
      "Analyte:", gene,
      "-", "Tumors:", paste(tumors, collapse = ", "),
      "- Stages:", paste(stages, collapse = ", ")
    )

    p <- ggplot(filtered_data, aes(x = .data$Tumor_Type, y = .data$z_score, fill = .data$Tumor_Type)) +
      geom_boxplot() +
      labs(title = title, x = "Tumor Type", y = "Abundance (Z-Score)") +
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

    if (isTRUE(show_points)) {
      p <- p +
        geom_point(
          aes(fill = .data$Tumor_Type),
          position = position_jitter(width = 0.2, height = 0),
          size = 1.5, alpha = 0.8,
          shape = 21, color = "#000000", stroke = 0.5
        )
    }

    if (isTRUE(show_significance)) {
      valid_tumors <- unique(filtered_data$Tumor_Type[filtered_data$Tumor_Type %in% tumors])
      unique_stages <- unique(filtered_data$Stage)

      # One tumor across multiple stages: compare stages
      if (length(valid_tumors) == 1 && length(unique_stages) > 1) {

        filtered_data$Tumor_Stage <- paste(filtered_data$Tumor_Type, filtered_data$Stage, sep = " - ")

        p <- ggplot(filtered_data, aes(x = .data$Tumor_Stage, y = .data$z_score, fill = .data$Tumor_Type)) +
          geom_boxplot() +
          labs(title = title, x = "Tumor Type - Stage", y = "Abundance (Z-Score)") +
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
          p <- p + stat_compare_means(
            method = "wilcox.test",
            comparisons = valid_comparisons,
            aes(group = Stage),
            label = "p.signif",
            size = 5
          )
        }

      } else if (length(valid_tumors) >= 2) {
        # Multiple tumors: compare tumors within each stage
        for (stage in unique_stages) {
          stage_data <- filtered_data[filtered_data$Stage == stage, ]
          stage_tumors <- unique(stage_data$Tumor_Type)

          if (length(stage_tumors) >= 2) {
            comparisons <- combn(stage_tumors, 2, simplify = FALSE)
            valid_comparisons <- comparisons[sapply(comparisons, function(x) {
              sum(stage_data$Tumor_Type == x[1]) > 0 & sum(stage_data$Tumor_Type == x[2]) > 0
            })]

            if (length(valid_comparisons) > 0) {
              p <- p + stat_compare_means(
                method = "wilcox.test",
                comparisons = valid_comparisons,
                aes(group = Tumor_Type),
                label = "p.signif",
                size = 5,
                data = stage_data
              )
            }
          }
        }
      }
    }

    return(p)
  }

  # ----------------------------
  # Tumor vs Non-tumor
  # ----------------------------
  if (comparison_type == "Tumor vs Non-tumor") {

    # Case: "all" stages -> collapse per Patient within each tumor
    if ("all" %in% stages) {

      filtered_data <- data %>%
        dplyr::filter(.data$Name == gene, .data$Tumor_Type %in% tumors) %>%
        dplyr::group_by(.data$Tumor_Type, .data$Patient_ID) %>%
        dplyr::summarise(z_score = mean(.data$z_score, na.rm = TRUE), .groups = "drop")

      if (nrow(filtered_data) == 0) {
        message("No data available for Tumor vs Non-tumor comparison.")
        return(NULL)
      }

      if (analyte_option == "significant") {
        filtered_data <- filtered_data %>% dplyr::filter(.data$z_score < -1.96 | .data$z_score > 1.96)
        if (nrow(filtered_data) == 0) {
          message("No significant analytes available for Tumor vs Non-tumor comparison.")
          return(NULL)
        }
      }

      filtered_data <- filtered_data %>%
        dplyr::mutate(Sample_Type = ifelse(grepl("\\.N$", .data$Patient_ID), "Non-tumor", .data$Tumor_Type))

      available_stages <- data %>%
        dplyr::filter(.data$Name == gene, .data$Tumor_Type %in% tumors, !grepl("\\.N$", .data$Patient_ID)) %>%
        dplyr::pull(.data$Stage) %>%
        unique() %>%
        na.omit() %>%
        sort()

      title <- paste0(
        "Analyte: ", gene,
        " – Tumor: ", paste(tumors, collapse = ", "),
        " – Stages: All (", paste(available_stages, collapse = ", "), ")"
      )

      p <- ggplot(filtered_data, aes(x = .data$Sample_Type, y = .data$z_score, fill = .data$Sample_Type)) +
        geom_boxplot() +
        labs(title = title, x = "Sample Type", y = "Abundance (Z-Score)") +
        scale_fill_manual(values = c("Non-tumor" = "gray", tumor_colors)) +
        theme_minimal() +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.line = element_line(color = "gray"),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 12, face = "bold")
        )

      if (isTRUE(show_points)) {
        p <- p +
          geom_point(
            aes(fill = .data$Sample_Type),
            position = position_jitter(width = 0.2, height = 0),
            size = 1.5, alpha = 0.8,
            shape = 21, color = "#000000", stroke = 0.5
          )
      }

      if (isTRUE(show_significance)) {
        valid_tumors <- unique(filtered_data$Tumor_Type)
        comparisons <- lapply(valid_tumors, function(tumor) c("Non-tumor", tumor))
        if (length(comparisons) > 0) {
          p <- p + stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif")
        }
      }

      return(p)
    }

    # Case: specific stages -> include Non-tumor samples + tumor samples in selected stages
    filtered_data <- data %>%
      dplyr::filter(.data$Name == gene, .data$Tumor_Type %in% tumors, .data$Stage %in% stages | grepl("\\.N$", .data$Patient_ID))

    if (nrow(filtered_data) == 0) {
      message("No data available for Tumor vs Non-tumor comparison.")
      return(NULL)
    }

    if (analyte_option == "significant") {
      filtered_data <- filtered_data %>% dplyr::filter(.data$z_score < -1.96 | .data$z_score > 1.96)
      if (nrow(filtered_data) == 0) {
        message("No significant analytes available for Tumor vs Non-tumor comparison.")
        return(NULL)
      }
    }

    filtered_data <- filtered_data %>%
      dplyr::mutate(Sample_Type = ifelse(grepl("\\.N$", .data$Patient_ID), "Non-tumor", paste(.data$Tumor_Type, .data$Stage, sep = " - ")))

    title <- paste(
      "Analyte:", gene,
      "- Tumor:", paste(tumors, collapse = ", "),
      "- Stages:", paste(stages, collapse = ", ")
    )

    p <- ggplot(filtered_data, aes(
      x = .data$Sample_Type,
      y = .data$z_score,
      fill = ifelse(.data$Sample_Type == "Non-tumor", "Non-tumor", .data$Tumor_Type)
    )) +
      geom_boxplot() +
      labs(title = title, x = "Sample Type", y = "Abundance (Z-Score)") +
      scale_fill_manual(values = c("Non-tumor" = "gray", tumor_colors)) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.line = element_line(color = "gray"),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 12, face = "bold")
      )

    if (isTRUE(show_points)) {
      p <- p +
        geom_point(
          aes(fill = .data$Sample_Type),
          position = position_jitter(width = 0.2, height = 0),
          size = 1.5, alpha = 0.8,
          shape = 21, color = "#000000", stroke = 0.5
        )
    }

    if (isTRUE(show_significance)) {
      valid_sample_types <- unique(filtered_data$Sample_Type)
      if (length(valid_sample_types) < 2) {
        message("Not enough valid sample types for significance testing.")
        return(NULL)
      }
      comparisons <- lapply(valid_sample_types[valid_sample_types != "Non-tumor"], function(sample) c("Non-tumor", sample))
      if (length(comparisons) > 0) {
        p <- p + stat_compare_means(method = "wilcox.test", comparisons = comparisons, label = "p.signif")
      } else {
        message("There are no valid tumor comparisons for significance testing.")
      }
    }

    return(p)
  }

  message("Unsupported comparison_type: ", comparison_type)
  return(NULL)
}

# =============================================================
# PROTEIN ACTIVITY: load choices and data
# =============================================================

load_protein_choices <- function() {
  gene_names_file <- "Prot_Act_names.tsv"
  if (!file.exists(gene_names_file)) {
    cat("Error: Gene names file not found: ", gene_names_file, "\n")
    return(NULL)
  }
  readr::read_tsv(gene_names_file, show_col_types = FALSE)$gene_name
}

load_data_protein <- function(protein) {
  gene_names_file <- "Prot_Act_names.tsv"
  if (!file.exists(gene_names_file)) {
    cat("Error: Protein names file not found: ", gene_names_file, "\n")
    return(NULL)
  }

  gene_names <- readr::read_tsv(gene_names_file, show_col_types = FALSE)$gene_name
  if (!protein %in% gene_names) {
    cat("Error: Protein not found in the gene names list.\n")
    return(NULL)
  }

  gene_letter <- substr(protein, 1, 1)
  file_name <- paste0("Protein_activity/Protein_activity_", gene_letter, ".tsv")
  cat("File path:", file_name, "\n")

  if (file.exists(file_name)) {
    return(readr::read_tsv(file_name, show_col_types = FALSE))
  }
  cat("Error: Protein Activity data file not found.\n")
  return(NULL)
}

# =============================================================
# PROTEIN ACTIVITY: plot generation
# =============================================================

generate_plot_protein_activity <- function(
    data,
    protein,
    tumors,
    stages,
    comparison_type_prot,
    show_points_protein = FALSE,
    show_significance = FALSE
) {
  if (is.null(data) || nrow(data) == 0) {
    message("No data provided to generate_plot_protein_activity().")
    return(NULL)
  }

  if (comparison_type_prot == "Tumor vs Tumor") {

    if ("all" %in% stages) {

      filtered_data <- data %>%
        dplyr::filter(.data$Name == protein, .data$Tumor_Type %in% tumors, !grepl("\\.N$", .data$Patient_ID)) %>%
        dplyr::group_by(.data$Tumor_Type, .data$Patient_ID) %>%
        dplyr::summarise(predicted_activity = mean(.data$predicted_activity, na.rm = TRUE), .groups = "drop")

      if (nrow(filtered_data) == 0) {
        message("No data available for the selected protein and tumors with 'all' stages.")
        return(NULL)
      }

      available_stages <- data %>%
        dplyr::filter(.data$Name == protein, .data$Tumor_Type %in% tumors, !grepl("\\.N$", .data$Patient_ID)) %>%
        dplyr::pull(.data$Stage) %>%
        unique() %>%
        na.omit() %>%
        sort()

      title <- paste0("Protein: ", protein,
                      " – Tumor: ", paste(tumors, collapse = ", "),
                      " – Stages: All (", paste(available_stages, collapse = ", "), ")")

      p <- ggplot(filtered_data, aes(x = .data$Tumor_Type, y = .data$predicted_activity, fill = .data$Tumor_Type)) +
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

      if (isTRUE(show_significance) && length(tumors) >= 2) {
        valid_tumors <- tumors[sapply(tumors, function(tumor) any(filtered_data$Tumor_Type == tumor))]
        if (length(valid_tumors) >= 2) {
          p <- p + stat_compare_means(method = "wilcox.test",
                                      comparisons = combn(valid_tumors, 2, simplify = FALSE),
                                      label = "p.signif")
        }
      }

      if (isTRUE(show_points_protein)) {
        p <- p + geom_point(aes(fill = .data$Tumor_Type),
                            position = position_jitter(width = 0.2, height = 0),
                            size = 1.5, alpha = 0.8,
                            shape = 21, color = "#000000", stroke = 0.5)
      }

      return(p)
    }

    filtered_data <- data %>%
      dplyr::filter(.data$Name == protein, .data$Tumor_Type %in% tumors, .data$Stage %in% stages, !grepl("\\.N$", .data$Patient_ID))

    if (nrow(filtered_data) == 0) {
      message("No data available for the selected protein, tumors, and stages.")
      return(NULL)
    }

    title <- paste("Protein:", protein, "-", "Tumors:", paste(tumors, collapse = ", "), "- Stages:", paste(stages, collapse = ", "))

    p <- ggplot(filtered_data, aes(x = .data$Tumor_Type, y = .data$predicted_activity, fill = .data$Tumor_Type)) +
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

    if (isTRUE(show_points_protein)) {
      p <- p + geom_point(aes(fill = .data$Tumor_Type),
                          position = position_jitter(width = 0.2, height = 0),
                          size = 1.5, alpha = 0.8,
                          shape = 21, color = "#000000", stroke = 0.5)
    }

    return(p)
  }

  if (comparison_type_prot == "Tumor vs Non-tumor") {

    if ("all" %in% stages) {
      filtered_data <- data %>%
        dplyr::filter(.data$Name == protein, .data$Tumor_Type %in% tumors) %>%
        dplyr::group_by(.data$Tumor_Type, .data$Patient_ID) %>%
        dplyr::summarise(predicted_activity = mean(.data$predicted_activity, na.rm = TRUE), .groups = "drop")

      if (nrow(filtered_data) == 0) {
        message("No data available for Tumor vs Non-tumor comparison.")
        return(NULL)
      }

      filtered_data <- filtered_data %>%
        dplyr::mutate(Sample_Type = ifelse(grepl("\\.N$", .data$Patient_ID), "Non-tumor", .data$Tumor_Type))

      title <- paste0("Protein: ", protein,
                      " – Tumor: ", paste(tumors, collapse = ", "),
                      " – Stages: All")

      p <- ggplot(filtered_data, aes(x = .data$Sample_Type, y = .data$predicted_activity, fill = .data$Sample_Type)) +
        geom_boxplot() +
        labs(title = title, x = "Sample Type", y = "Predicted activity") +
        scale_fill_manual(values = c("Non-tumor" = "gray", tumor_colors)) +
        theme_minimal() +
        theme(
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
          axis.text.y = element_text(size = 14),
          axis.line = element_line(color = "gray"),
          axis.title = element_text(size = 14),
          plot.title = element_text(size = 12, face = "bold")
        )

      if (isTRUE(show_points_protein)) {
        p <- p + geom_point(aes(fill = .data$Sample_Type),
                            position = position_jitter(width = 0.2, height = 0),
                            size = 1.5, alpha = 0.8,
                            shape = 21, color = "#000000", stroke = 0.5)
      }

      return(p)
    }

    filtered_data <- data %>%
      dplyr::filter(.data$Name == protein, .data$Tumor_Type %in% tumors, .data$Stage %in% stages | grepl("\\.N$", .data$Patient_ID))

    if (nrow(filtered_data) == 0) {
      message("No data available for Tumor vs Non-tumor comparison.")
      return(NULL)
    }

    filtered_data <- filtered_data %>%
      dplyr::mutate(Sample_Type = ifelse(grepl("\\.N$", .data$Patient_ID), "Non-tumor", paste(.data$Tumor_Type, .data$Stage, sep = " - ")))

    title <- paste("Analyte:", protein, "- Tumor:", paste(tumors, collapse = ", "), "- Stages:", paste(stages, collapse = ", "))

    p <- ggplot(filtered_data, aes(x = .data$Sample_Type,
                                  y = .data$predicted_activity,
                                  fill = ifelse(.data$Sample_Type == "Non-tumor", "Non-tumor", .data$Tumor_Type))) +
      geom_boxplot() +
      labs(title = title, x = "Sample Type", y = "Predicted activity") +
      scale_fill_manual(values = c("Non-tumor" = "gray", tumor_colors)) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.line = element_line(color = "gray"),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 12, face = "bold")
      )

    if (isTRUE(show_points_protein)) {
      p <- p + geom_point(aes(fill = .data$Sample_Type),
                          position = position_jitter(width = 0.2, height = 0),
                          size = 1.5, alpha = 0.8,
                          shape = 21, color = "#000000", stroke = 0.5)
    }

    return(p)
  }

  message("Unsupported comparison_type_prot: ", comparison_type_prot)
  return(NULL)
}

# =============================================================
# TABLE GENERATION (Tumor + matched Non-tumor)
# =============================================================

generate_table_norm_ab <- function(data, gene, tumor_type, stage, analyte_option = "all") {
  tumor_data <- data %>%
    dplyr::filter(.data$Name == gene) %>%
    dplyr::mutate(Type = dplyr::if_else(grepl("\\.N$", .data$Patient_ID), "Non-tumor", "Tumor"))

  tumor_data <- tumor_data %>%
    dplyr::filter(
      .data$Tumor_Type %in% as.character(tumor_type),
      .data$Type == "Non-tumor" | (.data$Type == "Tumor" & (is.na(.data$Stage) | "all" %in% stage | .data$Stage %in% stage))
    )

  if (identical(analyte_option, "significant")) {
    tumor_data <- tumor_data %>% dplyr::filter(.data$z_score < -1.96 | .data$z_score > 1.96)
  }

  tumor_data
}

generate_table_norm_ac <- function(data, protein, tumor_type, stage = "all") {
  tumor_data <- data %>%
    dplyr::filter(.data$Name == protein) %>%
    dplyr::mutate(Type = dplyr::if_else(grepl("\\.N$", .data$Patient_ID), "Non-tumor", "Tumor"))

  tumor_data <- tumor_data %>%
    dplyr::filter(
      .data$Tumor_Type %in% as.character(tumor_type),
      .data$Type == "Non-tumor" | (.data$Type == "Tumor" & (is.na(.data$Stage) | "all" %in% stage | .data$Stage %in% stage))
    )

  tumor_data
}

# =============================================================
# PTM network link (PhosphoSIGNOR)
# =============================================================

query_ptm_residue <- function(uniprot_id, protein_name, residue) {
  base_url <- "https://signor.uniroma2.it/PhosphoSIGNOR/apis/residueSearch.php"
  response <- httr::GET(base_url, query = list(id = uniprot_id, residue = residue))

  if (httr::http_error(response)) {
    cat("Error: Unable to connect to SIGNOR API.\n")
    return(NULL)
  }

  result <- httr::content(response, as = "text", encoding = "UTF-8")
  print(result)  # optional debug

  if (grepl("No result found!", result)) {
    cat("Residue not found in SIGNOR database.\n")
    return(NULL)
  }

  formatted_id <- paste0(uniprot_id, "_", protein_name, "_", residue)
  paste0("https://signor.uniroma2.it/PhosphoSIGNOR/results/entity.php?role=residue&ID=", formatted_id)
}

# =============================================================
# VIPER evidence helpers + PatientProfiler-style wrapper
# (kept as provided; comments are in English)
# =============================================================

.load_regulons_for_evidence <- function(analysis,
                                       organism,
                                       integrated_regulons = FALSE,
                                       collectri = FALSE,
                                       custom = FALSE,
                                       custom_path = NULL) {
  analysis <- match.arg(analysis, c("tfea", "ksea"))
  organism <- match.arg(organism, c("human", "mouse"))

  if (analysis == "ksea" && collectri) stop("Collectri is only a 'tfea' parameter.")
  if (custom && is.null(custom_path)) stop("Please provide a path to the custom regulon table.")

  if (custom) {
    message("Reading custom regulons (evidence)...")
    return(readr::read_tsv(custom_path, show_col_types = FALSE))
  }

  if (analysis == "tfea") {
    if (organism == "human") {
      if (collectri) {
        data("tfea_db_human_collectri", package = "SignalingProfiler", envir = environment())
        return(tfea_db_human_collectri)
      } else {
        data("tfea_db_human", package = "SignalingProfiler", envir = environment())
        return(tfea_db_human)
      }
    } else {
      data("tfea_db_mouse", package = "SignalingProfiler", envir = environment())
      return(tfea_db_mouse)
    }
  }

  if (organism == "human") {
    if (integrated_regulons) {
      return(access_remote_file(file = "ksea_db_human_atlas.tsv", dir = "PKN"))
    } else {
      data("ksea_db_human", package = "SignalingProfiler", envir = environment())
      return(ksea_db_human)
    }
  } else {
    data("ksea_db_mouse", package = "SignalingProfiler", envir = environment())
    return(ksea_db_mouse)
  }
}


build_viper_evidence_tables <- function(omic_data,
                                       analysis,
                                       organism,
                                       reg_minsize,
                                       exp_sign = FALSE,
                                       integrated_regulons = FALSE,
                                       collectri = FALSE,
                                       custom = FALSE,
                                       custom_path = NULL) {

  analysis <- match.arg(analysis, c("tfea", "ksea"))
  organism <- match.arg(organism, c("human", "mouse"))

  viper_format <- SignalingProfiler::create_viper_format(omic_data, analysis, significance = exp_sign)
  diff_matrix <- SignalingProfiler::create_matrix_from_VIPER_format(viper_format)

  regulons <- .load_regulons_for_evidence(
    analysis = analysis,
    organism = organism,
    integrated_regulons = integrated_regulons,
    collectri = collectri,
    custom = custom,
    custom_path = custom_path
  )

  if (!all(c("tf", "target") %in% colnames(regulons))) {
    stop("Regulon table must contain at least columns: 'tf' and 'target'.")
  }

  measured_targets <- rownames(diff_matrix)

  regulons_measured <- regulons %>%
    dplyr::filter(.data$target %in% measured_targets)

  if (nrow(regulons_measured) == 0) {
    stop("No measured analytes found in regulons. Check organism and target ID format (e.g., UNIPROT-aa-pos for KSEA).")
  }

  measured_by_tf <- regulons_measured %>%
    dplyr::group_by(.data$tf) %>%
    dplyr::summarise(
      n_targets_measured = dplyr::n_distinct(.data$target),
      targets_measured = list(sort(unique(.data$target))),
      .groups = "drop"
    )

  eligible_tfs <- measured_by_tf %>%
    dplyr::filter(.data$n_targets_measured >= reg_minsize) %>%
    dplyr::pull(.data$tf)

  regulons_used <- regulons_measured %>%
    dplyr::filter(.data$tf %in% eligible_tfs)

  regulons_used_by_tf <- measured_by_tf %>%
    dplyr::filter(.data$tf %in% eligible_tfs)

  if (nrow(regulons_used) == 0) {
    stop(paste0(
      "No regulators pass reg_minsize=", reg_minsize,
      " after intersecting regulons with measured targets."
    ))
  }

  regulons_total_by_tf <- regulons %>%
    dplyr::group_by(.data$tf) %>%
    dplyr::summarise(n_targets_total = dplyr::n_distinct(.data$target), .groups = "drop")

  if ("mor" %in% colnames(regulons_used)) {
    regulons_used <- regulons_used %>%
      dplyr::mutate(
        interaction_sign = dplyr::case_when(
          suppressWarnings(!is.na(as.numeric(.data$mor))) ~ as.numeric(.data$mor),
          tolower(as.character(.data$mor)) %in% c("inhibition", "inhibits", "repression", "repress") ~ -1,
          tolower(as.character(.data$mor)) %in% c("activation", "activates", "induction", "induces") ~ 1,
          TRUE ~ 1
        )
      )
  } else {
    regulons_used <- regulons_used %>% dplyr::mutate(interaction_sign = 1)
  }

  if ("weight" %in% colnames(regulons_used)) {
    regulons_used <- regulons_used %>%
      dplyr::mutate(interaction_weight = suppressWarnings(as.numeric(.data$weight)))
  } else {
    regulons_used <- regulons_used %>% dplyr::mutate(interaction_weight = 1)
  }
  regulons_used$interaction_weight[is.na(regulons_used$interaction_weight)] <- 1

  regulons_used <- regulons_used %>%
    dplyr::mutate(
      target_stat = as.numeric(diff_matrix[.data$target, 1]),
      interaction_score = interaction_sign * interaction_weight * target_stat
    ) %>%
    dplyr::left_join(regulons_total_by_tf, by = "tf") %>%
    dplyr::left_join(regulons_used_by_tf %>% dplyr::select(tf, n_targets_measured), by = "tf") %>%
    dplyr::mutate(
      coverage = dplyr::if_else(!is.na(n_targets_total) & n_targets_total > 0,
                                n_targets_measured / n_targets_total,
                                NA_real_)
    )

  return(list(
    regulons_used = regulons_used,
    regulons_used_by_tf = regulons_used_by_tf
  ))
}


run_footprint_based_analysis_evidence <- function(omic_data,
                                                 analysis,
                                                 organism,
                                                 reg_minsize,
                                                 exp_sign = FALSE,
                                                 integrated_regulons = FALSE,
                                                 collectri = FALSE,
                                                 hypergeom_corr = TRUE,
                                                 correct_proteomics = FALSE,
                                                 prot_df = NULL,
                                                 GO_annotation = FALSE,
                                                 custom = FALSE,
                                                 custom_path = NULL) {

  out <- run_footprint_based_analysis(
    omic_data = omic_data,
    analysis = analysis,
    organism = organism,
    reg_minsize = reg_minsize,
    exp_sign = exp_sign,
    integrated_regulons = integrated_regulons,
    collectri = collectri,
    hypergeom_corr = hypergeom_corr,
    correct_proteomics = correct_proteomics,
    prot_df = prot_df,
    GO_annotation = GO_annotation,
    custom = custom,
    custom_path = custom_path
  )

  ev <- build_viper_evidence_tables(
    omic_data = omic_data,
    analysis = analysis,
    organism = organism,
    reg_minsize = reg_minsize,
    exp_sign = exp_sign,
    integrated_regulons = integrated_regulons,
    collectri = collectri,
    custom = custom,
    custom_path = custom_path
  )

  attr(out, "regulons_used") <- ev$regulons_used
  attr(out, "regulons_used_by_tf") <- ev$regulons_used_by_tf
  return(out)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

extract_patient_id <- function(filename) {
  x <- gsub("\\.(tsv|csv|xlsx|xls)$", "", filename, ignore.case = TRUE)
  x <- gsub("^Prot_Patient_|^Transc_Patient_|^Phospho_Patient_", "", x)
  x <- gsub("^P_|^T_|^Ph_", "", x)
  x
}

extract_protein_activity <- function(
    prot_file = NULL,
    trans_file = NULL,
    phospho_file = NULL,
    tf_params = list(),
    kin_params = list(),
    phosphoscore_params = list(),
    phosphoscore_noseqwin_params = list(),
    output_dir = "Activities",
    return_regulons = TRUE,
    save_regulons = FALSE,
    regulons_dir = file.path(output_dir, "Regulons")
) {

  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (isTRUE(save_regulons) && !dir.exists(regulons_dir)) dir.create(regulons_dir, recursive = TRUE)

  Prot_P <- if (!is.null(prot_file)) readr::read_tsv(prot_file, show_col_types = FALSE) else NULL
  Trans_P <- if (!is.null(trans_file)) readr::read_tsv(trans_file, show_col_types = FALSE) else NULL
  Phospho_P <- if (!is.null(phospho_file)) readr::read_tsv(phospho_file, show_col_types = FALSE) else NULL

  if (!is.null(Phospho_P) && "UNIPROT" %in% colnames(Phospho_P)) {
    Phospho_P <- Phospho_P %>%
      dplyr::mutate(UNIPROT = as.character(.data$UNIPROT)) %>%
      tidyr::separate_rows(.data$UNIPROT, sep = ";") %>%
      dplyr::mutate(UNIPROT = trimws(.data$UNIPROT)) %>%
      dplyr::filter(!is.na(.data$UNIPROT), .data$UNIPROT != "")

    key_cols <- intersect(c("UNIPROT", "aminoacid", "position", "gene_name", "difference", "logpval"), colnames(Phospho_P))
    if (length(key_cols) >= 2) {
      Phospho_P <- Phospho_P %>% dplyr::distinct(dplyr::across(dplyr::all_of(key_cols)), .keep_all = TRUE)
    }
  }

  default_tf_params <- list(
    omic_data = Trans_P,
    analysis = "tfea",
    organism = "human",
    reg_minsize = 10,
    exp_sign = FALSE,
    collectri = FALSE,
    hypergeom_corr = TRUE,
    GO_annotation = TRUE,
    correct_proteomics = FALSE,
    prot_df = Prot_P,
    custom = FALSE,
    custom_path = NULL
  )
  tf_params <- modifyList(default_tf_params, tf_params)

  default_kin_params <- list(
    omic_data = Phospho_P,
    analysis = "ksea",
    organism = "human",
    reg_minsize = 5,
    exp_sign = FALSE,
    integrated_regulons = TRUE,
    hypergeom_corr = TRUE,
    GO_annotation = TRUE,
    correct_proteomics = FALSE,
    prot_df = Prot_P,
    custom = FALSE,
    custom_path = NULL
  )
  kin_params <- modifyList(default_kin_params, kin_params)

  default_phosphoscore_params <- list(
    phosphoproteomic_data = Phospho_P,
    organism = "human",
    activatory = TRUE,
    GO_annotation = TRUE,
    custom = FALSE,
    custom_path = NULL
  )
  phosphoscore_params <- modifyList(default_phosphoscore_params, phosphoscore_params)

  default_phosphoscore_noseqwin_params <- list(
    phosphoproteomic_data = Phospho_P,
    organism = "human",
    activatory = TRUE,
    GO_annotation = TRUE,
    custom = FALSE,
    custom_path = NULL
  )
  phosphoscore_noseqwin_params <- modifyList(default_phosphoscore_noseqwin_params, phosphoscore_noseqwin_params)

  tf_activity_foot_1 <- if (!is.null(Trans_P)) {
    do.call(SignalingProfiler::run_footprint_based_analysis, tf_params)
  } else data.frame()

  kin_phos_activity_foot_1 <- if (!is.null(Phospho_P)) {
    do.call(SignalingProfiler::run_footprint_based_analysis, kin_params)
  } else data.frame()

  regulons_out <- list(
    tf = list(regulons_used = NULL, regulons_used_by_tf = NULL),
    kin = list(regulons_used = NULL, regulons_used_by_tf = NULL)
  )

  if (isTRUE(return_regulons) && exists("run_footprint_based_analysis_evidence", mode = "function")) {
    if (!is.null(Trans_P)) {
      tf_ev <- tryCatch(do.call(run_footprint_based_analysis_evidence, tf_params), error = function(e) NULL)
      if (!is.null(tf_ev)) {
        regulons_out$tf$regulons_used <- attr(tf_ev, "regulons_used")
        regulons_out$tf$regulons_used_by_tf <- attr(tf_ev, "regulons_used_by_tf")
      }
    }
    if (!is.null(Phospho_P)) {
      kin_ev <- tryCatch(do.call(run_footprint_based_analysis_evidence, kin_params), error = function(e) NULL)
      if (!is.null(kin_ev)) {
        regulons_out$kin$regulons_used <- attr(kin_ev, "regulons_used")
        regulons_out$kin$regulons_used_by_tf <- attr(kin_ev, "regulons_used_by_tf")
      }
    }
  }

  toy_activity_df <- data.frame()

  if (!is.null(Phospho_P)) {

    if (!("sequence_window" %in% colnames(Phospho_P))) {
      Phospho_P$sequence_window <- NA_character_
    }

    if (!all(is.na(Phospho_P$sequence_window))) {
      phosphoscore_1 <- do.call(SignalingProfiler::phosphoscore_computation, phosphoscore_params)
      other_method <- "PhosphoScore"
    } else {
      phosphoscore_1 <- do.call(SignalingProfiler::phosphoscore_computation_aapos, phosphoscore_noseqwin_params)
      other_method <- "PhosphoScore_AAPOS"
    }

    combined_tf <- if (nrow(tf_activity_foot_1) > 0) {
      SignalingProfiler::combine_footprint_and_phosphoscore(
        footprint_output = tf_activity_foot_1,
        phosphoscore_df = phosphoscore_1,
        analysis = "tfea"
      )
    } else data.frame()

    combined_kin_phos <- if (nrow(kin_phos_activity_foot_1) > 0) {
      SignalingProfiler::combine_footprint_and_phosphoscore(
        footprint_output = kin_phos_activity_foot_1,
        phosphoscore_df = phosphoscore_1,
        analysis = "ksea"
      )
    } else data.frame()

    toy_other <- phosphoscore_1 %>%
      dplyr::filter(.data$mf == "other") %>%
      dplyr::rename(final_score = .data$phosphoscore) %>%
      dplyr::mutate(method = other_method)

    toy_phos <- phosphoscore_1 %>%
      dplyr::filter(.data$mf == "phos") %>%
      dplyr::rename(final_score = .data$phosphoscore) %>%
      dplyr::mutate(method = "PhosphoScore")

    toy_activity_df <- dplyr::bind_rows(combined_tf, combined_kin_phos, toy_phos, toy_other) %>%
      dplyr::select(.data$UNIPROT, .data$gene_name, .data$mf, .data$final_score, .data$method)

  } else {

    toy_activity_df <- dplyr::bind_rows(tf_activity_foot_1, kin_phos_activity_foot_1) %>%
      dplyr::mutate(method = "VIPER") %>%
      dplyr::rename(final_score = .data$weightedNES) %>%
      dplyr::select(.data$UNIPROT, .data$gene_name, .data$mf, .data$final_score, .data$method)
  }

  if (!is.null(toy_activity_df) && nrow(toy_activity_df) > 0) {
    patient_name <- extract_patient_id(basename(prot_file %||% trans_file %||% phospho_file))
    readr::write_tsv(toy_activity_df, file.path(output_dir, paste0("Activity_Patient_", patient_name, ".tsv")))

    if (isTRUE(save_regulons) && isTRUE(return_regulons)) {
      if (is.data.frame(regulons_out$tf$regulons_used)) {
        readr::write_tsv(regulons_out$tf$regulons_used, file.path(regulons_dir, paste0("Regulons_TF_used_", patient_name, ".tsv")))
      }
      if (is.data.frame(regulons_out$kin$regulons_used)) {
        readr::write_tsv(regulons_out$kin$regulons_used, file.path(regulons_dir, paste0("Regulons_KIN_used_", patient_name, ".tsv")))
      }
    }
  } else {
    warning("No activity data could be extracted.")
  }

  if (isTRUE(return_regulons)) {
    return(list(activity = toy_activity_df, regulons = regulons_out))
  } else {
    return(invisible(toy_activity_df))
  }
}
