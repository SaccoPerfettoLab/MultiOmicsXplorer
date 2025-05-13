library(SignalingProfiler)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(dplyr)

#### READ UPLOADED FILES ####

read_data <- function(filepath) {
  ext <- tools::file_ext(filepath)
  
  switch(ext,
         csv = {
          
           first_line <- readLines(filepath, n = 1) # read first line to guess delimiter -> cause can be csv or csv2!!
           
           if (grepl(";", first_line)) {
             readr::read_csv2(filepath)  # semicolon-separated
           } else {
             readr::read_csv(filepath)   # comma-separated
           }
         },
         tsv = readr::read_tsv(filepath),
         xlsx = readxl::read_xlsx(filepath),
         xls = readxl::read_xls(filepath),
         stop("File type not supported.")
  )
}

#### Manage NA character problem ####
# suggested by vero to manage if the na is character
normalize_na <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.character(col)) {
      col[col == "NA"] <- NA
    }
    return(col)
  })
  return(df)
}

#### RUN ANALYSIS ####
perform_analysis <- function(trans_data = NULL, prot_data = NULL, phospho_data = NULL, params, sample_id) {

  organism <- params$organism
  
  Tras_P <- if (!is.null(trans_data) && !is.data.frame(trans_data)) normalize_na(read_data(trans_data)) else trans_data
  Prot_P <- if (!is.null(prot_data) && !is.data.frame(prot_data)) normalize_na(read_data(prot_data)) else prot_data
  Phospho_P <- if (!is.null(phospho_data) && !is.data.frame(phospho_data)) normalize_na(read_data(phospho_data)) else phospho_data
  
  if (is.null(Tras_P) && is.null(Phospho_P)) {
    warning("No transcriptomic or phosphoproteomic data uploaded for patient: ", sample_id)
    return(list(combined_results = NULL, molecular_function_summary = NULL))
  }
  
  toy_activity_df <- data.frame()
  
  tf_activity_foot_1 <- NULL
  kin_phos_activity_foot_1 <- NULL
  combined_results <- list()  
  
  if (!is.null(Tras_P)) {
    tf_args <- list(
      omic_data = Tras_P, 
      analysis = 'tfea', 
      organism = organism,
      reg_minsize = params$reg_minsize_tf, 
      exp_sign = params$exp_sign_tf,
      hypergeom_corr = params$hypergeom_corr_tf, 
      GO_annotation = params$GO_annotation_tf,
      collectri = FALSE, 
      correct_proteomics = params$correct_proteomics_tf
    )
    
    if (params$correct_proteomics_tf) {
      if (is.null(Prot_P)) {
        stop("Proteomics data is required for TFEA correction, but none provided.")
      }
      tf_args$prot_df <- Prot_P  
    }
    
    tf_activity_foot_1 <- do.call(SignalingProfiler::run_footprint_based_analysis, tf_args)
    
    if (!is.null(tf_activity_foot_1)) {
      tf_activity_foot_1 <- tf_activity_foot_1 %>%
        mutate(method = if (!"method" %in% colnames(.)) "VIPER" else method)
      combined_results <- append(combined_results, list(tf_activity_foot_1))
    }
  }
  
  if (!is.null(Phospho_P)) {
     ## partial cleaning
    Phospho_P <- Phospho_P %>%
      filter(aminoacid %in% c("S", "T", "Y")) %>% ## if there are aminoacids that are not T, Y ans S, the rows containing them will be deleted
      distinct(UNIPROT, aminoacid, position, .keep_all = TRUE) # it mantains only one row per uniprot-aminoacid_position combination -> this could be not correct if a gene_name 
                                                              # as the wrong uniprot cause it mantains the first occurrance 
    
    kin_args <- list(
      omic_data = Phospho_P, 
      analysis = 'ksea', 
      organism = organism,
      reg_minsize = params$reg_minsize_kin, 
      exp_sign = params$exp_sign_kin,
      hypergeom_corr = params$hypergeom_corr_kin, 
      GO_annotation = params$GO_annotation_kin,
      correct_proteomics = params$correct_proteomics_kin
    )
    
    if (params$correct_proteomics_kin) {
      if (is.null(Prot_P)) {
        message("Proteomics data is required for KSEA correction, but none provided.")
      }
      kin_args$prot_df <- Prot_P  
    }
    
    kin_phos_activity_foot_1 <- do.call(SignalingProfiler::run_footprint_based_analysis, kin_args)
    
    
    if (organism == "human") {  
      phosphoscore_1 <- if (!all(is.null(Phospho_P$sequence_window))) {
        SignalingProfiler::phosphoscore_computation(
          phosphoproteomic_data = Phospho_P, organism = 'human',
          activatory = params$activatory, GO_annotation = params$GO_annotation_phospho
        )
      } else {
        SignalingProfiler::phosphoscore_computation_aapos(
          phosphoproteomic_data = Phospho_P, organism = 'human',
          activatory = params$activatory, GO_annotation = params$GO_annotation_phospho
        )
      }
    } else {  
      phosphoscore_1 <- if (!all(is.null(Phospho_P$sequence_window))) {
        SignalingProfiler::phosphoscore_computation(
          phosphoproteomic_data = Phospho_P, organism = 'hybrid', ## to allineate the sequence on human sequences
          activatory = params$activatory, GO_annotation = params$GO_annotation_phospho,blastp_path = "/Users/eleonorameo/bin/ncbi-blast-2.16.0+/bin/blastp" #re-write if blastp path change
        )
      } else {
        SignalingProfiler::phosphoscore_computation_aapos(
          phosphoproteomic_data = Phospho_P, organism = 'hybrid', 
          activatory = params$activatory, GO_annotation = params$GO_annotation_phospho
        )
      }
    }
    
    
    if ("mf" %in% colnames(phosphoscore_1)) {
      toy_other <- phosphoscore_1 %>%
        dplyr::filter(mf == 'other') %>%
        dplyr::rename(final_score = phosphoscore) %>%
        dplyr::mutate(method = "PhosphoScore")
    } else {
      toy_other <- phosphoscore_1 %>%
        dplyr::rename(final_score = phosphoscore) %>%
        dplyr::mutate(method = "PhosphoScore")
    }
    
    if (!is.null(tf_activity_foot_1)) {
      combined_tf <- SignalingProfiler::combine_footprint_and_phosphoscore(
        footprint_output = tf_activity_foot_1, phosphoscore_df = phosphoscore_1, analysis = 'tfea'
      )
      combined_results <- append(combined_results, list(combined_tf))
    }
    
    if (!is.null(kin_phos_activity_foot_1)) {
      combined_kin_phos <- SignalingProfiler::combine_footprint_and_phosphoscore(
        footprint_output = kin_phos_activity_foot_1, phosphoscore_df = phosphoscore_1, analysis = 'ksea'
      )
      combined_results <- append(combined_results, list(combined_kin_phos))
    }
    
    combined_results <- append(combined_results, list(toy_other))  
  }
  
  if (length(combined_results) > 0) {
    final_combined_df <- bind_rows(combined_results) %>%
      mutate(UNIPROT = sub("\\;.*", "", UNIPROT)) %>%
      unique()
    
    if ("final_score" %in% colnames(final_combined_df)) {
      final_results <- final_combined_df %>%
        select(any_of(c("UNIPROT", "gene_name", "mf", "final_score", "method"))) %>%
        rename(predicted_activity = final_score)
    } else if ("weightedNES" %in% colnames(final_combined_df)) {
      final_results <- final_combined_df %>%
        select(any_of(c("UNIPROT", "gene_name", "mf", "weightedNES", "method"))) %>%
        rename(predicted_activity = weightedNES)
    } else {
      final_results <- final_combined_df %>%
        select(any_of(c("UNIPROT", "gene_name", "mf"))) %>%
        mutate(predicted_activity = NA, method = NA)
    }
    
    
    final_results <- final_results %>%
      filter(!is.na(predicted_activity) & predicted_activity != 0) %>%
      mutate(Sample_ID = sample_id)
  } else {
    final_results <- data.frame(UNIPROT = character(), gene_name = character(),
                                mf = character(), predicted_activity = numeric(), 
                                method = character(), Sample_ID = character())
  }
  
  ## create the table with the molecular function summary
  molecular_function_summary <- if (nrow(final_results) > 0 && "mf" %in% colnames(final_results)) {
    final_results %>%
      group_by(Sample_ID) %>%
      summarise(
        tf_count = sum(mf == "tf", na.rm = TRUE),
        kin_count = sum(mf == "kin", na.rm = TRUE),
        phos_count = sum(mf == "phos", na.rm = TRUE),
        other_count = sum(mf == "other", na.rm = TRUE)
      )
  } else {
    data.frame(Sample_ID = sample_id, tf_count = 0, kin_count = 0, phos_count = 0, other_count = 0)
  }
  
  return(list(
    combined_results = final_results,
    molecular_function_summary = molecular_function_summary
  ))
}



#### PLOT RESULTS ####
create_top_proteins_plot <- function(results_df, sample_id = NULL) {

  if (!is.null(sample_id)) {
    cat("Filtering for Sample ID:", sample_id, "\n")
    results_df <- results_df %>% filter(Sample_ID == !!sample_id)  
    cat("Rows after filtering:", nrow(results_df), "\n")
  }
  
  top_proteins <- data.frame()
  
  categories <- unique(results_df$mf)
  
  for (cat in categories) {
    filtered_df <- results_df %>% filter(mf == cat)  
    
    if (nrow(filtered_df) > 0) {
      if (is.null(sample_id)) {
        filtered_df <- filtered_df %>%
          group_by(gene_name, mf) %>%
          summarize(predicted_activity = mean(predicted_activity, na.rm = TRUE), .groups = 'drop')
      }
      
      top_up <- filtered_df %>%
        filter(predicted_activity > 0) %>%  
        arrange(desc(predicted_activity)) %>%
        head(15) %>%  
        mutate(Type = "Up-regulated")
      
      top_down <- filtered_df %>%
        filter(predicted_activity < 0) %>%
        arrange(predicted_activity) %>%
        head(15) %>%  
        mutate(Type = "Down-regulated")
      
      top_cat <- bind_rows(top_up, top_down)
      
      top_proteins <- bind_rows(top_proteins, top_cat)
    }
  }
  
  
  # re-writed to be similar to the plot of signaling profiler tutorial
  # vertical plot
  plot <- ggplot(top_proteins, aes(x = reorder(gene_name, predicted_activity), 
                                   y = predicted_activity, 
                                   fill = predicted_activity)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3, linetype = "solid") + 
    labs(
      title = if (is.null(sample_id)) {
        "Top 15 Up and Down-regulated Proteins for each molecular function (All Samples)"
      } else {
        paste("Top 15 Up and Down-regulated Proteins for each molecular function - Sample ID:", sample_id)
      },
      x = "Protein Name", 
      y = "Predicted Activity", 
      fill = "Predicted Activity"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.position = "bottom",
      strip.text = element_text(size = 12, face = "bold")
    ) +
    scale_fill_gradient2(
      low = "#4F9DFF", mid = "white", high = "#FF6F61",
      midpoint = 0,
      limits = c(-10, 10),
      oob = scales::squish,
      breaks = c(-10, 0, 10),
      name = "Predicted Activity"
    ) +
    scale_y_discrete(expand = expansion(mult = c(0.05, 0.05))) +
    facet_wrap(~ mf, scales = "free_y") +
    coord_flip()
  
  
  
  return(plot)
}



# #### CHECK ERRORS ####
check_columns <- function(data, omic_type) {
  if (omic_type == "Transcriptomics") {
    required_cols <- c("gene_name", "difference", "logpval", "significant")
  } else if (omic_type == "Proteomics") {
    required_cols <- c("gene_name", "UNIPROT", "difference", "logpval", "significant")
  } else if (omic_type == "Phosphoproteomics") {
    required_cols <- c("UNIPROT", "gene_name", "aminoacid", "position", "difference", "logpval", "significant")
  } else {
    return(FALSE)
  }

  missing_cols <- setdiff(required_cols, colnames(data)) # check if there's not the right columns

  if (length(missing_cols) > 0) {
    return(paste("The input table must include at least the following columns for", omic_type, "data: ",
                 paste(required_cols, collapse = ", ")))
  }

  return(TRUE)
}

