
# ---------- Optional helper for Shiny UI choices ----------
if (!exists("load_protein_choices")) {
  load_protein_choices <- function() {
    gene_names_file <- "Prot_Act_names.tsv"
    if (!file.exists(gene_names_file)) {
      cat("Error: Gene names file not found:", gene_names_file, "\n")
      return(NULL)
    }
    readr::read_tsv(gene_names_file, show_col_types = FALSE)$gene_name
  }
}

# ---------- Libraries ----------
suppressPackageStartupMessages({
  library(SignalingProfiler)
  library(dplyr)
  library(tidyverse)
  library(ggplot2)
})

# ============================================================
# 1) FILE READING UTILITIES
# ============================================================

# Read uploaded file by extension and auto-detect CSV delimiter.
read_data <- function(filepath) {
  ext <- tools::file_ext(filepath)

  switch(ext,
         csv = {
           # Read the first line to guess the delimiter (',' vs ';')
           first_line <- readLines(filepath, n = 1)

           if (grepl(";", first_line)) {
             readr::read_csv2(filepath)  # semicolon-separated
           } else {
             readr::read_csv(filepath)   # comma-separated
           }
         },
         tsv  = readr::read_tsv(filepath, show_col_types = FALSE),
         xlsx = readxl::read_xlsx(filepath),
         xls  = readxl::read_xls(filepath),
         stop("File type not supported: ", ext)
  )
}

# Convert string "NA" to actual NA (common issue in uploads).
normalize_na <- function(df) {
  df[] <- lapply(df, function(col) {
    if (is.character(col)) {
      col[col == "NA"] <- NA
    }
    col
  })
  df
}

# ============================================================
# 2) INPUT TABLE VALIDATION (Shiny uploader requirement)
# ============================================================

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

  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    return(paste(
      "The input table must include at least the following columns for",
      omic_type, "data:",
      paste(required_cols, collapse = ", ")
    ))
  }
  TRUE
}

# ============================================================
# 3) MAIN ANALYSIS FUNCTION
# ============================================================

perform_analysis <- function(trans_data = NULL,
                             prot_data = NULL,
                             phospho_data = NULL,
                             params,
                             sample_id) {

  organism <- params$organism

  # ---- Load / normalize inputs ----
  Tras_P   <- if (!is.null(trans_data) && !is.data.frame(trans_data)) normalize_na(read_data(trans_data)) else trans_data
  Prot_P   <- if (!is.null(prot_data)  && !is.data.frame(prot_data))  normalize_na(read_data(prot_data))  else prot_data
  Phospho_P <- if (!is.null(phospho_data) && !is.data.frame(phospho_data)) normalize_na(read_data(phospho_data)) else phospho_data

  # Nothing to do
  if (is.null(Tras_P) && is.null(Phospho_P)) {
    warning("No transcriptomic or phosphoproteomic data uploaded for patient: ", sample_id)
    return(list(
      combined_results = NULL,
      molecular_function_summary = NULL,
      regulons = data.frame()
    ))
  }

  # ---- Initialize containers ----
  tf_activity_foot_1       <- NULL
  kin_phos_activity_foot_1 <- NULL
  phosphoscore_1           <- NULL
  combined_results         <- list()

  # Evidence containers (regulons)
  regulons_tables  <- list()
  regulons_summary <- list()

  # ==========================================================
  # 3A) TFEA (Transcriptomics) - VIPER footprint
  # ==========================================================
  if (!is.null(Tras_P)) {

    tf_args <- list(
      omic_data           = Tras_P,
      analysis            = "tfea",
      organism            = organism,
      reg_minsize         = params$reg_minsize_tf,
      exp_sign            = params$exp_sign_tf,
      hypergeom_corr      = params$hypergeom_corr_tf,
      GO_annotation       = params$GO_annotation_tf,
      collectri           = FALSE,
      correct_proteomics  = params$correct_proteomics_tf
    )

    if (isTRUE(params$correct_proteomics_tf)) {
      if (is.null(Prot_P)) stop("Proteomics data is required for TFEA correction, but none provided.")
      tf_args$prot_df <- Prot_P
    }

    tf_activity_foot_1 <- do.call(SignalingProfiler::run_footprint_based_analysis, tf_args)

    # Ensure a method label exists
    if (!is.null(tf_activity_foot_1)) {
      tf_activity_foot_1 <- tf_activity_foot_1 %>%
        dplyr::mutate(method = if (!"method" %in% colnames(.)) "VIPER" else method)
    }

    # ---- Evidence extraction (does NOT affect scores) ----
    tf_reg_used <- NULL
    tf_reg_sum  <- NULL

    if (exists("run_footprint_based_analysis_evidence", mode = "function")) {
      tf_args_ev <- tf_args
      tf_args_ev$integrated_regulons <- FALSE
      tf_ev <- tryCatch(do.call(run_footprint_based_analysis_evidence, tf_args_ev),
                        error = function(e) NULL)
      if (!is.null(tf_ev)) {
        tf_reg_used <- attr(tf_ev, "regulons_used")
        tf_reg_sum  <- attr(tf_ev, "regulons_used_by_tf")
      }
    } else {
      # Fallback: some versions may attach attributes to the main output
      tf_reg_used <- attr(tf_activity_foot_1, "regulons_used")
      tf_reg_sum  <- attr(tf_activity_foot_1, "regulons_used_by_tf")
    }

    if (!is.null(tf_reg_used) && nrow(tf_reg_used) > 0) {
      tf_reg_used <- tf_reg_used %>%
        dplyr::mutate(Sample_ID = sample_id, analysis = "tfea")
      regulons_tables <- append(regulons_tables, list(tf_reg_used))
    }

    if (!is.null(tf_reg_sum) && nrow(tf_reg_sum) > 0) {
      tf_reg_sum <- tf_reg_sum %>%
        dplyr::mutate(Sample_ID = sample_id, analysis = "tfea")
      regulons_summary <- append(regulons_summary, list(tf_reg_sum))
    }
  }

  # If only transcriptomics is provided, return TFEA-only results
  if (!is.null(tf_activity_foot_1) && is.null(Phospho_P)) {
    combined_results <- append(combined_results, list(tf_activity_foot_1))
  }

  # ==========================================================
  # 3B) KSEA (Phosphoproteomics footprint) + PhosphoScore
  # ==========================================================
  if (!is.null(Phospho_P)) {

    # Basic cleaning: keep only S/T/Y sites
    Phospho_P <- Phospho_P %>% dplyr::filter(aminoacid %in% c("S", "T", "Y"))

    kin_args <- list(
      omic_data           = Phospho_P,
      analysis            = "ksea",
      organism            = organism,
      reg_minsize         = params$reg_minsize_kin,
      exp_sign            = params$exp_sign_kin,
      integrated_regulons = params$integrated_regulons_kin,
      hypergeom_corr      = params$hypergeom_corr_kin,
      GO_annotation       = params$GO_annotation_kin,
      correct_proteomics  = params$correct_proteomics_kin
    )

    if (isTRUE(params$correct_proteomics_kin)) {
      if (is.null(Prot_P)) {
        message("Proteomics data is required for KSEA correction, but none provided. Disabling correction.")
        kin_args$correct_proteomics <- FALSE
      } else {
        kin_args$prot_df <- Prot_P
      }
    }

    kin_phos_activity_foot_1 <- do.call(SignalingProfiler::run_footprint_based_analysis, kin_args)

    # Evidence extraction (does NOT affect scores)
    kin_reg_used <- NULL
    kin_reg_sum  <- NULL

    if (exists("run_footprint_based_analysis_evidence", mode = "function")) {
      kin_args_ev <- kin_args
      kin_ev <- tryCatch(do.call(run_footprint_based_analysis_evidence, kin_args_ev),
                         error = function(e) NULL)
      if (!is.null(kin_ev)) {
        kin_reg_used <- attr(kin_ev, "regulons_used")
        kin_reg_sum  <- attr(kin_ev, "regulons_used_by_tf")
      }
    } else {
      kin_reg_used <- attr(kin_phos_activity_foot_1, "regulons_used")
      kin_reg_sum  <- attr(kin_phos_activity_foot_1, "regulons_used_by_tf")
    }

    if (!is.null(kin_reg_used) && nrow(kin_reg_used) > 0) {
      kin_reg_used <- kin_reg_used %>%
        dplyr::mutate(Sample_ID = sample_id, analysis = "ksea")
      regulons_tables <- append(regulons_tables, list(kin_reg_used))
    }

    if (!is.null(kin_reg_sum) && nrow(kin_reg_sum) > 0) {
      kin_reg_sum <- kin_reg_sum %>%
        dplyr::mutate(Sample_ID = sample_id, analysis = "ksea")
      regulons_summary <- append(regulons_summary, list(kin_reg_sum))
    }

    # Ensure a method label exists
    if (!is.null(kin_phos_activity_foot_1)) {
      kin_phos_activity_foot_1 <- kin_phos_activity_foot_1 %>%
        dplyr::mutate(method = if (!"method" %in% colnames(.)) "VIPER" else method)
    }

    # ----------------------------------------------------------
    # PhosphoScore computation (optional, fail-safe)
    # ----------------------------------------------------------
    phosphoscore_1 <- tryCatch({
      if (organism == "human") {
        if ("sequence_window" %in% colnames(Phospho_P) && !all(is.na(Phospho_P$sequence_window))) {
          SignalingProfiler::phosphoscore_computation(
            phosphoproteomic_data = Phospho_P, organism = "human",
            activatory = params$activatory, GO_annotation = params$GO_annotation_phospho
          )
        } else {
          SignalingProfiler::phosphoscore_computation_aapos(
            phosphoproteomic_data = Phospho_P, organism = "human",
            activatory = params$activatory, GO_annotation = params$GO_annotation_phospho
          )
        }
      } else {
        if ("sequence_window" %in% colnames(Phospho_P) && !all(is.na(Phospho_P$sequence_window))) {
          SignalingProfiler::phosphoscore_computation(
            phosphoproteomic_data = Phospho_P, organism = "hybrid",
            activatory = params$activatory, GO_annotation = params$GO_annotation_phospho
          )
        } else {
          SignalingProfiler::phosphoscore_computation_aapos(
            phosphoproteomic_data = Phospho_P, organism = "hybrid",
            activatory = params$activatory, GO_annotation = params$GO_annotation_phospho
          )
        }
      }
    }, error = function(e) {
      message("PhosphoScore failed for patient ", sample_id, ": ", conditionMessage(e))
      NULL
    })

    # If PhosphoScore is missing/empty, return footprint-only results
    if (is.null(phosphoscore_1) || nrow(phosphoscore_1) == 0) {

      if (!is.null(tf_activity_foot_1)) {
        combined_results <- append(combined_results, list(tf_activity_foot_1))
      }
      if (!is.null(kin_phos_activity_foot_1)) {
        combined_results <- append(combined_results, list(kin_phos_activity_foot_1))
      }

    } else {

      # PhosphoScore outputs:
      # - mf == "phos"  : phosphatases activity
      # - mf == "other" : other protein classes
      if ("mf" %in% colnames(phosphoscore_1)) {
        toy_other <- phosphoscore_1 %>%
          dplyr::filter(.data$mf == "other") %>%
          dplyr::rename(final_score = .data$phosphoscore) %>%
          dplyr::mutate(method = "PhosphoScore")

        toy_phos <- phosphoscore_1 %>%
          dplyr::filter(.data$mf == "phos") %>%
          dplyr::rename(final_score = .data$phosphoscore) %>%
          dplyr::mutate(method = "PhosphoScore")
      } else {
        # If mf is missing, keep all PhosphoScore results as "other"
        toy_other <- phosphoscore_1 %>%
          dplyr::rename(final_score = .data$phosphoscore) %>%
          dplyr::mutate(mf = "other", method = "PhosphoScore")
        toy_phos <- NULL
      }

      # Combine TF footprint + phosphoscore
      if (!is.null(tf_activity_foot_1)) {
        combined_tf <- SignalingProfiler::combine_footprint_and_phosphoscore(
          footprint_output = tf_activity_foot_1,
          phosphoscore_df   = phosphoscore_1,
          analysis          = "tfea"
        )
        if (is.null(combined_tf) || nrow(combined_tf) == 0) combined_tf <- tf_activity_foot_1
        combined_results <- append(combined_results, list(combined_tf))
      }

      # Combine KIN footprint + phosphoscore
      if (!is.null(kin_phos_activity_foot_1)) {
        combined_kin_phos <- SignalingProfiler::combine_footprint_and_phosphoscore(
          footprint_output = kin_phos_activity_foot_1,
          phosphoscore_df   = phosphoscore_1,
          analysis          = "ksea"
        )
        if (is.null(combined_kin_phos) || nrow(combined_kin_phos) == 0) combined_kin_phos <- kin_phos_activity_foot_1
        combined_results <- append(combined_results, list(combined_kin_phos))
      }

      # Add phosphoscore-only classes
      combined_results <- append(combined_results, list(toy_other))
      if (!is.null(toy_phos) && nrow(toy_phos) > 0) combined_results <- append(combined_results, list(toy_phos))
    }
  }

  # ==========================================================
  # 3C) FINAL COMBINE (PatientProfiler-compatible)
  # ==========================================================

  # Drop NULL elements (bind_rows would fail)
  combined_results <- combined_results[!vapply(combined_results, is.null, logical(1))]

  if (length(combined_results) > 0) {

    final_combined_df <- dplyr::bind_rows(combined_results)

    # Standardize the activity column name to `predicted_activity`
    if ("final_score" %in% colnames(final_combined_df)) {
      final_results <- final_combined_df %>%
        dplyr::select(dplyr::any_of(c("UNIPROT", "gene_name", "mf", "final_score", "method"))) %>%
        dplyr::rename(predicted_activity = .data$final_score)
    } else if ("weightedNES" %in% colnames(final_combined_df)) {
      final_results <- final_combined_df %>%
        dplyr::select(dplyr::any_of(c("UNIPROT", "gene_name", "mf", "weightedNES", "method"))) %>%
        dplyr::rename(predicted_activity = .data$weightedNES)
    } else if ("NES" %in% colnames(final_combined_df)) {
      final_results <- final_combined_df %>%
        dplyr::select(dplyr::any_of(c("UNIPROT", "gene_name", "mf", "NES", "method"))) %>%
        dplyr::rename(predicted_activity = .data$NES)
    } else {
      final_results <- final_combined_df %>%
        dplyr::select(dplyr::any_of(c("UNIPROT", "gene_name", "mf", "method"))) %>%
        dplyr::mutate(predicted_activity = NA_real_)
    }

    # IMPORTANT: do NOT filter-out NA/0 here; PatientProfiler keeps the full set
    final_results <- final_results %>% dplyr::mutate(Sample_ID = sample_id)

  } else {
    final_results <- data.frame(
      UNIPROT = character(),
      gene_name = character(),
      mf = character(),
      predicted_activity = numeric(),
      method = character(),
      Sample_ID = character()
    )
  }

  molecular_function_summary <- data.frame(
    Sample_ID = sample_id,
    tf_count = sum(final_results$mf == "tf", na.rm = TRUE),
    kin_count = sum(final_results$mf == "kin", na.rm = TRUE),
    phos_count = sum(final_results$mf == "phos", na.rm = TRUE),
    other_count = sum(final_results$mf == "other", na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  # ==========================================================
  # 3D) BUILD EVIDENCE TABLE: "regulons"
  # ==========================================================

  final_regulons_used <- dplyr::bind_rows(regulons_tables)
  if (is.null(final_regulons_used) || nrow(final_regulons_used) == 0) {
    final_regulons_used <- data.frame()
  }

  # Filter regulons to keep only those that plausibly contributed to inferred activities:
  #  - keep only regulators that appear in inferred activities (final_results)
  #  - keep only targets that are measurable in the provided omics tables
  #  - if target-level stat exists, drop NA (targets not used because not in signature)
  if (nrow(final_regulons_used) > 0 && !is.null(final_results) && nrow(final_results) > 0) {

    inferred_tf  <- unique(final_results$gene_name[final_results$mf == "tf"])
    inferred_kin <- unique(final_results$gene_name[final_results$mf == "kin"])

    measurable_tfea <- if (!is.null(Tras_P) && "gene_name" %in% colnames(Tras_P)) {
      unique(as.character(Tras_P$gene_name))
    } else character()

    measurable_ksea <- if (!is.null(Phospho_P) && all(c("UNIPROT", "aminoacid", "position") %in% colnames(Phospho_P))) {
      Phospho_P %>%
        tidyr::separate_rows(UNIPROT, sep = ";") %>%
        dplyr::mutate(target_id = paste0(UNIPROT, "-", aminoacid, "-", position)) %>%
        dplyr::pull(.data$target_id) %>%
        unique()
    } else character()

    final_regulons_used <- final_regulons_used %>%
      dplyr::mutate(
        analysis      = if ("analysis" %in% colnames(.)) as.character(.data$analysis) else NA_character_,
        regulator_id  = if ("tf" %in% colnames(.)) as.character(.data$tf) else NA_character_,
        target_id     = if ("target" %in% colnames(.)) as.character(.data$target) else NA_character_
      ) %>%
      dplyr::filter(
        dplyr::case_when(
          .data$analysis == "tfea" ~ .data$regulator_id %in% inferred_tf,
          .data$analysis == "ksea" ~ .data$regulator_id %in% inferred_kin,
          TRUE ~ FALSE
        )
      ) %>%
      dplyr::filter(
        dplyr::case_when(
          .data$analysis == "tfea" ~ .data$target_id %in% measurable_tfea,
          .data$analysis == "ksea" ~ .data$target_id %in% measurable_ksea,
          TRUE ~ TRUE
        )
      )

    if ("target_stat" %in% colnames(final_regulons_used)) {
      final_regulons_used <- final_regulons_used %>% dplyr::filter(!is.na(.data$target_stat))
    } else if ("stat" %in% colnames(final_regulons_used)) {
      final_regulons_used <- final_regulons_used %>% dplyr::filter(!is.na(.data$stat))
    }
  }

  # Build KSEA target mapping from the input phosphoproteomics table:
  # UNIPROT-AA-POS  ->  GENE_AApos (e.g., O15511-S-77 -> ARPC5_S77)
  ksea_target_map <- NULL
  if (!is.null(Phospho_P) && all(c("UNIPROT", "aminoacid", "position", "gene_name") %in% colnames(Phospho_P))) {
    ksea_target_map <- Phospho_P %>%
      tidyr::separate_rows(UNIPROT, sep = ";") %>%
      dplyr::mutate(
        target_id = paste0(UNIPROT, "-", aminoacid, "-", position),
        target_gene_site = paste0(gene_name, "_", aminoacid, position)
      ) %>%
      dplyr::select(target_id, target_gene_site) %>%
      dplyr::distinct()
  }

  regulons_viper <- if (nrow(final_regulons_used) > 0) {

    interaction_sign <- if ("interaction_sign" %in% colnames(final_regulons_used)) {
      suppressWarnings(as.numeric(final_regulons_used$interaction_sign))
    } else if ("mor" %in% colnames(final_regulons_used)) {
      suppressWarnings(as.numeric(final_regulons_used$mor))
    } else {
      rep(1, nrow(final_regulons_used))
    }
    interaction_sign[is.na(interaction_sign)] <- 1

    stat <- if ("target_stat" %in% colnames(final_regulons_used)) {
      suppressWarnings(as.numeric(final_regulons_used$target_stat))
    } else if ("stat" %in% colnames(final_regulons_used)) {
      suppressWarnings(as.numeric(final_regulons_used$stat))
    } else {
      rep(NA_real_, nrow(final_regulons_used))
    }

    tmp <- final_regulons_used %>%
      dplyr::mutate(
        sample_id = .data$Sample_ID,
        `protein estimated in activity` = .data$tf,
        target_id = .data$target,
        `analysis type` = if ("analysis" %in% colnames(.)) as.character(.data$analysis) else NA_character_,
        `source of information type` = "regulon",
        `regulatory effect of source of information` = interaction_sign,
        `stat of source information` = stat
      )

    # Map KSEA target IDs to gene-site labels if possible
    if (!is.null(ksea_target_map)) {
      tmp <- tmp %>%
        dplyr::left_join(ksea_target_map, by = c("target_id" = "target_id")) %>%
        dplyr::mutate(
          `source of information` = dplyr::if_else(
            !is.na(.data$`analysis type`) &
              .data$`analysis type` == "ksea" &
              !is.na(.data$target_gene_site),
            .data$target_gene_site,
            .data$target_id
          )
        )
    } else {
      tmp <- tmp %>% dplyr::mutate(`source of information` = .data$target_id)
    }

    tmp %>%
      dplyr::mutate(
        `overall effect on protein activity` =
          suppressWarnings(as.numeric(.data$`regulatory effect of source of information`)) *
          suppressWarnings(as.numeric(.data$`stat of source information`))
      ) %>%
      dplyr::select(
        sample_id,
        `protein estimated in activity`,
        `source of information`,
        `source of information type`,
        `analysis type`,
        `regulatory effect of source of information`,
        `stat of source information`,
        `overall effect on protein activity`
      ) %>%
      dplyr::distinct()

  } else {
    data.frame(
      sample_id = character(),
      `protein estimated in activity` = character(),
      `source of information` = character(),
      `source of information type` = character(),
      `analysis type` = character(),
      `regulatory effect of source of information` = numeric(),
      `stat of source information` = numeric(),
      `overall effect on protein activity` = numeric(),
      check.names = FALSE
    )
  }

  regulons_phosphoscore <- if (!is.null(phosphoscore_1) && nrow(phosphoscore_1) > 0) {

    has_phos <- "phos" %in% colnames(phosphoscore_1)
    has_act  <- "act_rol" %in% colnames(phosphoscore_1)
    has_fc   <- "phosphosite_fc" %in% colnames(phosphoscore_1)

    if (!(has_phos && has_act && has_fc)) {
      data.frame(
        sample_id = character(),
        `protein estimated in activity` = character(),
        `source of information` = character(),
        `source of information type` = character(),
        `analysis type` = character(),
        `regulatory effect of source of information` = numeric(),
        `stat of source information` = numeric(),
        `overall effect on protein activity` = numeric(),
        check.names = FALSE
      )
    } else {

      phosphoscore_1 %>%
        tidyr::separate_rows(phos, act_rol, phosphosite_fc, sep = ";") %>%
        dplyr::mutate(
          phos = stringr::str_trim(phos),
          act_rol = stringr::str_trim(act_rol),
          phosphosite_fc = stringr::str_trim(phosphosite_fc),
          sample_id = sample_id,
          `protein estimated in activity` = .data$gene_name,
          `source of information` = paste0(.data$gene_name, "_", phos),
          `source of information type` = "phosphosite",
          `analysis type` = "phosphoscore",
          `regulatory effect of source of information` = suppressWarnings(as.numeric(.data$act_rol)),
          `stat of source information` = suppressWarnings(as.numeric(.data$phosphosite_fc)),
          `overall effect on protein activity` =
            suppressWarnings(as.numeric(.data$act_rol)) *
            suppressWarnings(as.numeric(.data$phosphosite_fc))
        ) %>%
        dplyr::select(
          sample_id,
          `protein estimated in activity`,
          `source of information`,
          `source of information type`,
          `analysis type`,
          `regulatory effect of source of information`,
          `stat of source information`,
          `overall effect on protein activity`
        ) %>%
        dplyr::distinct()
    }

  } else {
    data.frame(
      sample_id = character(),
      `protein estimated in activity` = character(),
      `source of information` = character(),
      `source of information type` = character(),
      `analysis type` = character(),
      `regulatory effect of source of information` = numeric(),
      `stat of source information` = numeric(),
      `overall effect on protein activity` = numeric(),
      check.names = FALSE
    )
  }

  regulons <- dplyr::bind_rows(regulons_viper, regulons_phosphoscore) %>% dplyr::distinct()

  return(list(
    combined_results = final_results,
    molecular_function_summary = molecular_function_summary,
    regulons = regulons
  ))
}

# ============================================================
# 4) PLOTTING (Top up/down per molecular function)
# ============================================================

create_top_proteins_plot <- function(results_df, sample_id = NULL, top_n = 15) {

  if (!is.null(sample_id)) {
    cat("Filtering for Sample ID:", sample_id, "\n")
    results_df <- results_df %>% dplyr::filter(Sample_ID == !!sample_id)
    cat("Rows after filtering:", nrow(results_df), "\n")
  }

  top_n_label <- if (is.null(top_n) || is.infinite(top_n)) "All" else as.character(top_n)

  top_proteins <- data.frame()
  categories <- unique(results_df$mf)

  for (cat in categories) {
    filtered_df <- results_df %>% dplyr::filter(mf == cat)

    if (nrow(filtered_df) > 0) {

      # If not filtering by a specific sample, aggregate by mean across samples
      if (is.null(sample_id)) {
        filtered_df <- filtered_df %>%
          dplyr::group_by(gene_name, mf) %>%
          dplyr::summarize(predicted_activity = mean(predicted_activity, na.rm = TRUE), .groups = "drop")
      }

      top_up <- filtered_df %>%
        dplyr::filter(predicted_activity > 0) %>%
        dplyr::arrange(desc(predicted_activity)) %>%
        head(top_n) %>%
        dplyr::mutate(Type = "Up-regulated")

      top_down <- filtered_df %>%
        dplyr::filter(predicted_activity < 0) %>%
        dplyr::arrange(predicted_activity) %>%
        head(top_n) %>%
        dplyr::mutate(Type = "Down-regulated")

      top_cat <- dplyr::bind_rows(top_up, top_down)
      top_proteins <- dplyr::bind_rows(top_proteins, top_cat)
    }
  }

  ggplot(top_proteins, aes(x = reorder(gene_name, predicted_activity),
                           y = predicted_activity,
                           fill = predicted_activity)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3, linetype = "solid") +
    labs(
      title = if (is.null(sample_id)) {
        paste0("Top ", top_n_label, " Up and Down-regulated Proteins for each molecular function (All Samples)")
      } else {
        paste0("Top ", top_n_label, " Up and Down-regulated Proteins for each molecular function - Sample ID: ", sample_id)
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
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    facet_wrap(~ mf, scales = "free_y") +
    coord_flip()
}
