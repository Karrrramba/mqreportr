# Functions

# Load packages
# library(BiocManager)
# library(missForest) #amputation
# library(randomForest) #imputation
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
  print(quant_regex)
  d <- readr::read_tsv(file = data,
                       col_select = c(
                         cols,
                         matches(quant_regex)
                         )
                       ) %>% 
    janitor::clean_names(., abbreviation = "ID")
  
  # Remove non-normalized intensity/ratio columns if applicable
  d<- d[, !(grepl("^(reporter|lfq|ratio)", colnames(d)) & !grepl("(corrected|normalized)", colnames(d)))]

  return(d)
}

### Universal data cleaning----
# Removes reverse, contaminants and on;y id by site, filters sites to desired threshold
clean_df <- function(data, prob_treshold = 0.9) {
  if("only" %in%  names(data)){
    oisr <- which(data[, "only_identified_by_site"] == "+")
    oisc <- which("only_identified_by_site" %in% names(data))
    data <- data[-oisr, -oisc]
  } 
  # Remove decoy entries and conatminants
  del_row <- which(data[, "reverse"] == "+" | data[, "potential_contaminant"] == "+")
  del_col <- which(names(data) %in% c("reverse", "potential_contaminant"))
  data <- data[-del_row, -del_col]
  data <- data["localization_prob" >= prob_treshold, ]
  
  # Collapse protein ids and gene names to first entry
  if("protein_id_s" %in% names(data)){
    data <- rename(data, "protein" = "protein_id_s")
  }
  
  data <- data %>% 
    mutate(
      protein = gsub("(;.+)", "", protein),
      gene_names = gsub("(;.+)", "", gene_names)
    )
  
  return(data)
}

# Transform into long format----
# Helper function
remove_empty_channels <- function(data){
  ch <- data %>% 
    group_by(channel) %>% 
    summarise(values = length(unique(log2_intensity))) %>% 
    filter(values <= 2) %>% 
    pull(channel)
  
  data <- filter(data, !channel %in% ch)
  
  return(data)
}

transform_table <- function(data, site = FALSE) {
  data <- data %>% 
    tidyr::pivot_longer(
      cols = matches("^(reporter|lfq|ratio|intensity)"), 
      names_to = "sample",
      values_to = "log2_intensity"
    ) %>% 
    dplyr::mutate(
      # Remove prefix
      sample = trimws(gsub("(reporter|lfq|ratio|intensity|corrected)", "", sample), which = "left", whitespace = "_"),
      channel = as.numeric(stringr::str_extract(sample, "^\\d+")),
      sample = gsub("^(\\d+_)", "", sample)
    ) %>% 
    # Log2-transform
    dplyr::mutate(log2_intensity = dplyr::if_else(log2_intensity != 0, log2(log2_intensity), 0))
  
  if ("channel" %in% names(data)){
    data <- remove_empty_channels(data)
  }
  
  if (site == TRUE){
    data <- data %>% 
      mutate(
        multiplicity = stringr::str_extract(sample, "\\d$"),
        sample = gsub("_1-3]", "", sample),
        site = paste0(amino_acid, position)
      ) %>% 
      select(!c(amino_acid, position))
  }
  return(data)
}


# Experimental Design Templates----
# Creates an annotation file in the working directory
write_annotation_file <- function(data) {
  if ("channel" %in% names(data)){
    a <- data %>% 
      dplyr::group_by(sample) %>% 
      dplyr::reframe(channel = unique(channel)) %>% 
      dplyr::mutate(group1 = NA)
  } else {
  a <- data.frame(sample = unique(data$sample),
                  group_1 = NA)
  }
  
  write.csv(a, file = "annotation_file.csv", row.names = FALSE)
}

# Read annotation file and create groups based on column names
load_annotations <- function(annotation, data) {
  annotations <- readr::read_csv(a_file)
  
  # Merge annotations with data to create groups
  data <- data %>% 
    dplyr::left_join(., annotations, by = c("sample", "channel"), relationship = "many-to-many")
  
  return(data)
}

