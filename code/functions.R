# Functions

# Load packages
# library(BiocManager)
# library(missForest) #amputation
# library(randomForest) #imputation
library(data.table)
library(dtplyr)
library(tidyverse)

## Parse files----
# mq_tables <- list.files(mq_output, pattern = '\\.txt$')
# site_tables <- NULL
# for (t in mq_tables) {
#   if (str_detect(t, "Sites")){
#     site_tables <- append(site_tables, t)
#   }
# }

# Load all detected files into global env
# list2env(
#   lapply(setNames(site_tables, janitor::make_clean_names(site_tables)), load_table, multi = "lfq"),
#   envir = .GlobalEnv)

# Table loading----
remove_prefix <- function(data){
  names(data) <- trimws(gsub("(reporter|lfq|ratio|intensity|corrected)", "", names(data)), which = "left", whitespace = "_")
  return(data)
}

load_table <- function(data, quant) {
  
  cols <- c(
    "Protein names", 
    "Gene names", 
    "Reverse",
    "Potential contaminant"
    )
  
  pg_cols <- c(
    "Protein IDs",
    "Only identified by site"
    )
  
  site_cols <- c(
    "Protein",
    "Amino acid", 
    "Position",
    "Localization prob",
    "Sequence window"
    )
  
  # Add modification-relevant columns if applicable 
  if (stringr::str_detect(data, "Sites")) {
    cols <- c(cols, site_cols)
  } else if(stringr::str_detect(data, "Groups")) {
    cols <- c(cols, pg_cols)
  }

  # Select quantification columns based on quantification strategy
  if (quant == "lfq") {
    quant_regex <- "^(LFQ )|(I|i)ntensity .+"
  } else if (quant == "tmt") {
    quant_regex <- "Reporter intensity"
  } else {
    quant_regex <- "Ratio"
  }

  d <- dtplyr::lazy_dt(readr::read_tsv(file = data,
                       col_select = c(
                         cols,
                         tidyselect::matches(quant_regex)
                         )
                       ) %>% 
    janitor::clean_names(., abbreviation = "ID")) %>% 
    data.table::as.data.table()
  
  # Remove non-normalized intensity/ratio columns if applicable
  # d <- d[, !(grepl("^(reporter|lfq|ratio)", colnames(d)) & !grepl("(corrected|normalized)", colnames(d)))]
  d <- d %>% dtplyr::lazy_dt() %>% 
    dplyr::select(!(tidyselect::starts_with("Reporter") & !tidyselect::matches("corrected"))) %>% 
    data.table::as.data.table()

  return(d)
}

### Universal data cleaning----
# Removes reverse, contaminants and on;y id by site, filters sites to desired threshold
clean_df <- function(data, prob_treshold = 0.9) {
  if("only" %in%  names(data)){
    data <- data %>% 
      dtplyr::lazy_dt() %>% 
      dplyr::filter(is.na(only_identified_by_site)) %>% 
      dplyr::select(!only_identified_by_site) %>% 
      data.table::as.data.table()
  } 
  # Remove decoy entries and contaminants
  data <- data %>% 
    dtplyr::lazy_dt() %>%
    dplyr::filter(is.na(reverse) & is.na(potential_contaminant)) %>% 
    dplyr::filter(localization_prob >= prob_treshold) %>% 
    dplyr::select(!c(reverse, potential_contaminant)) %>% 
    dplyr::mutate(
      protein = gsub("(;.+)", "", protein),
      gene_names = gsub("(;.+)", "", gene_names)
    ) %>% 
    data.table::as.data.table()
  
  # Collapse protein ids and gene names to first entry
  if("protein_id_s" %in% names(data)){
    data <- dplyr::rename(data, "protein" = "protein_id_s")
  }
  
  return(data)
}

# Transform into long format----
# Helper function
remove_empty_channels <- function(data){
  ch <- data %>% 
    dtplyr::lazy_dt() %>% 
    dplyr::group_by(channel) %>% 
    dplyr::summarise(values = length(dplyr::distinct(log2_intensity))) %>% 
    dplyr::filter(values <= 2) %>% 
    dplyr::pull(channel) %>% 
    data.table::as.data.table()
  
  data <- dplyr::filter(data, !channel %in% ch) #needs adjustment either to base R or dtplyr
  
  return(data)
}

transform_long <- function(data, site = FALSE) {
  data <- data %>% 
    dtplyr::lazy_dt() %>% 
    tidyr::pivot_longer(
      cols = tidyselect::matches("^(reporter|lfq|ratio|intensity)"), 
      names_to = "sample",
      values_to = "log2_intensity"
    ) %>% 
    dplyr::mutate(
      # Remove prefix
      sample = trimws(gsub("(reporter|lfq|ratio|intensity|corrected)", "", sample), which = "left", whitespace = "_"),
      # Extract channel from sample name (not applicable for non-tmt as-is)
      channel = as.numeric(stringr::str_extract(sample, "^\\d+")),
      sample = gsub("^(\\d+_)", "", sample)
    ) %>% 
    # Log2-transformation 
    dplyr::mutate(log2_intensity = dplyr::if_else(log2_intensity != 0, log2(log2_intensity), NA)) %>% 
    dplyr::relocate(channel, .after = sample) %>% 
    data.table::as.data.table()
  
  if ("channel" %in% names(data)){
    data <- remove_empty_channels(data)
  }
  # Add sites and multiplicity for PTMs 
  if (site == TRUE){
    data <- data %>% 
      dtplyr::lazy_dt() %>% 
      dplyr::mutate(
        multiplicity = stringr::str_extract(sample, "\\d$"),
        sample = gsub("_[1-3]$", "", sample),
        site = paste0(amino_acid, position)
      ) %>% 
      dplyr::select(!c(amino_acid, position)) %>% 
      dplyr::relocate(c(site, multiplicity), .after = protein) %>% 
      data.table::as.data.table()
  }
  return(data)
}

# Experimental Design Templates----
# Creates an annotation file in the working directory
write_group_annotations_file <- function(data) {
  if ("channel" %in% names(data)){
    out <- data %>% 
      dplyr::group_by(sample) %>% 
      dplyr::reframe(channel = unique(channel)) %>% 
      dplyr::mutate(group1 = NA)
  } else {
  out <- data.frame(sample = unique(data$sample),
                  group_1 = NA)
  }
  
  write.csv(
    out, 
    file = paste0(
      deparse(substitute(data)), 
      "_annotations_", 
      format(Sys.time(), "%Y-%m-%d %H.%M"), 
      ".csv"
      ), 
    row.names = FALSE
  )
}

# Read annotation file and create groups based on column names
load_group_annotations <- function(annotation_file, data) {
  annotations <- readr::read_csv(annotation_file)
  
  # Merge annotations with data to create groups
  data <- data %>% 
    dtplyr::lazy_dt() %>% 
    dplyr::left_join(., annotations, by = c("sample", "channel"), relationship = "many-to-many") %>% 
    data.table::as.data.table()
  
  return(data)
}

