# Functions

# Load packages
# library(BiocManager)
# library(COMBAT)
# library(limma) #loess normalization
# library(missForest) #amputation
# library(mixomics) #NIPALS PCA
# library(randomForest) #imputation
library(tidyverse)

## Parse files
mq_output <- "data-raw/" #user input
mq_tables <- list.files(mq_output, pattern = '\\.txt$')
site_tables <- NULL
for (t in mq_tables) {
  if (str_detect(t, "Site")){
    site_tables <- append(site_tables, t)
  }
}

remove_prefix <- function(data){
  names(data) <- trimws(gsub("(reporter|lfq|ratio|intensity|corrected)", "", names(data)), which = "left", whitespace = "_")
  return(data)
}

load_table <- function(data, quant) {
  
  cols <- c("Protein IDs",
            "Protein names", 
            "Gene names", 
            "Reverse",
            "Potential contaminant",
            "Only identified by site")
  
  site_cols <- c("Amino acid", 
                 "Site",
                 "Localization prob",
                 "Multiplicity",
                 "Sequence window")
  
  # Add modification-relevant columns if applicable 
  if (stringr::str_detect(data, "site")) {
    cols <- c(cols, site_cols)
  }
  
  # Select quantification columns based on quantification strategy
  if (quant == "lfq") {
    quant_regex <- "^(LFQ) |(I|i)ntensity .*"
  } else if (quant == "tmt") {
    quant_regex <- "Reporter intensity (corrected|\\d)"
  } else {
    quant_regex <- "Ratio .*"
  }
  
  d <- readr::read_tsv(file = data,
                        col_select = c(
                          cols,
                          matches(quant_regex))
                        ) %>% 
    janitor::clean_names(., abbreviation = "ID")
  
  # Remove non-normalized intensity/ratio columns if applicable
  d_clean <- d[, !(grepl("^(reporter|lfq|ratio)", colnames(d)) & !grepl("(corrected|normalized)", colnames(d)))]

  return(d_clean)
}

# Load all detected files into global env
list2env(
  lapply(setNames(site_tables, janitor::make_clean_names(site_tables)), load_table, multi = "lfq"),
  envir = .GlobalEnv)

### Universal data cleaning----
# Removes reverse, contaminants and on;y id by site, filters sites to desired threshold
clean_df <- function(data, prob_treshold = 0.9, remove_empty = TRUE) {
  del_row <- which(data[, "reverse"] == "+" | data[, "potential_contaminant"] == "+" | data[, "only_identified_by_site"] == "+")
  del_col <- which(names(data) %in% c("reverse", "potential_contaminant", "only_identified_by_site"))
  
  data <- data[-del_row, -del_col]
  data <- data["localization_prob" >= prob_treshold, ]
  
  if (remove_empty) {
    data <- remove_empty_channels(data)
  }
  
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

transform_table <- function(data) {
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
      sample = gsub("^(\\d+_)", "", sample),
      # Collapse ID columns to first entry
      protein_id_s = gsub("(;.+)", "", protein_id_s),
      gene_names = gsub("(;.+)", "", gene_names),
    ) %>% 
    # Log2-transform
    dplyr::mutate(log2_intensity = dplyr::if_else(log2_intensity != 0, log2(log2_intensity), 0))
  
  if ("channel" %in% names(data)){
    data <- remove_empty_channels(data)
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

rat_pg_annotated <- load_annotations(a_file = "annotation_file.csv", rat_pg_long) 


# proteinGroups----
rat_pg <- load_table('data-raw/tmt_rat/proteinGroups.txt', quant = "tmt")
rat_pg <- clean_df(rat_pg)
rat_pg_long <- transform_table(rat_pg)

rat_pg_nc <- rat_pg_long %>% 
  filter(grepl("nc", sample, fixed = TRUE))

write_annotation_file(rat_pg_nc)

rat_pg_annotated <- load_annotations(annotation = "annotation_file.csv", data = rat_pg_nc)


# sites



nem_sites <- readr::read_tsv("NEM (C)Sites.txt",
                             col_select =  c(
                               "Protein",
                               "Gene names",
                               "Protein names",
                               matches("Intensity .*"),
                               "Localization prob",
                               "Sequence window",
                               "Position",
                               "Reverse",
                               "Potential contaminant"
                             ))
nem_sites <- clean_data(nem_sites)
