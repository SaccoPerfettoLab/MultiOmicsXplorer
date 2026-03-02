
############################################################
#  Install required packages 
############################################################

# Install SignalingProfiler
devtools::install_github('https://github.com/SaccoPerfettoLab/SignalingProfiler/')

# Install PatientProfiler
devtools::install_github('https://github.com/SaccoPerfettoLab/PatientProfiler/')




############################################################
# Load required libraries
############################################################

library(PatientProfiler)
library(SignalingProfiler)
library(tidyverse)


############################################################
# 2) Download and load example datasets
############################################################

# Load phosphoproteomics data
url <- "https://perfettolab.bio.uniroma1.it/PerfettoLabData/PatientProfiler/MOX_example/Phosphoproteomics_updated.tsv"
file_path <- tempfile(fileext = ".tsv")
response <- httr::GET(url, httr::write_disk(file_path, overwrite = TRUE))
if (httr::http_status(response)$category == "Success") {
  ph_updated <- utils::read.delim(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
} else {
  stop("Error:", httr::http_status(response)$message)
}

# Load proteomics data
url <- "https://perfettolab.bio.uniroma1.it/PerfettoLabData/PatientProfiler/MOX_example/Proteomics_updated.tsv"
file_path <- tempfile(fileext = ".tsv")
response <- httr::GET(url, httr::write_disk(file_path, overwrite = TRUE))
if (httr::http_status(response)$category == "Success") {
  pr_updated <- utils::read.delim(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
} else {
  stop("Error:", httr::http_status(response)$message)
}

# Load trasncriptomics data
url <- "https://perfettolab.bio.uniroma1.it/PerfettoLabData/PatientProfiler/MOX_example/Transcriptomics_updated.tsv"
file_path <- tempfile(fileext = ".tsv")
response <- httr::GET(url, httr::write_disk(file_path, overwrite = TRUE))
if (httr::http_status(response)$category == "Success") {
  tr_updated <- utils::read.delim(file_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
} else {
  stop("Error:", httr::http_status(response)$message)
}


############################################################
# 3) Data preparation for patient-specific activity inference analysis
############################################################

# The function omics_preparation() converts the multi-patient
# updated omics tables into separate patient-specific files.
#
# For each sample, it:
# - extracts the molecular measurements (stored as "difference")
# - computes a p-value
# - assigns a significance flag based on a z-score threshold (±1.96)

# For additional informations about this function, you can check the 
# PatientProfiler GitHub https://github.com/SaccoPerfettoLab/PatientProfiler


omics_preparation(
  df_tr_updated = tr_updated,
  df_pr_updated = pr_updated,
  df_ph_updated = ph_updated,
  transc_dir_name = "Mox_data",
  prot_dir_name = "Mox_data",
  phospho_dir_name = "Mox_data"
)

# This function generates patient-specific files such as:
#
#   Phospho_Patient_CPT000814.tsv
#   Prot_Patient_CPT000814.tsv
#   Transc_Patient_CPT000814.tsv


# Rename files to Extract protein activity from your data/ Run analysis

data_dir <- "Mox_data"
sample_id <- "CPT000814"

file.rename(
  from = file.path(data_dir, "Prot_Patient_CPT000814.tsv"),
  to   = file.path(data_dir, paste0("P_", sample_id, ".tsv"))
)

file.rename(
  from = file.path(data_dir, "Phospho_Patient_CPT000814.tsv"),
  to   = file.path(data_dir, paste0("Ph_", sample_id, ".tsv"))
)

file.rename(
  from = file.path(data_dir, "Transc_Patient_CPT000814.tsv"),
  to   = file.path(data_dir, paste0("T_", sample_id, ".tsv"))
)

















