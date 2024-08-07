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

# Removes reverse and contaminants, filters sites to desired threshold
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

# Creates an annotation file in the working directory
write_annotation_file <- function(data) {
  
  nq_cols <- c("protein_id_s",
            "protein_names", 
            "gene_names", 
            "amino_acid", 
            "site",
            "localization_prob",
            "multiplicity",
            "sequence_window")
  
  a <- data.frame(sample = names(data[, !names(data) %in% nq_cols]),
                  group_1 = NA)
  
  write.csv(a, file = "annotation_file.csv")
}


# Read annotation file and create groups based on column names
read_annotation_file <- function(a_file, data) {
  
  annotations <- readr::read_csv(a_file)
  
  # Remove channels if applicable
  annotations <- annotations["keep" == TRUE, ]
  names(data) <- annotations[, names(data) %in% annotations$sample] 
  
  # Merge annotations with data to create groups
  data <- data %>% 
    dplyr::left_join()
  
  return(data)
}





# Experimental Design Templates----





# proteinGroups----
rat_pg <- load_table('data-raw/tmt_rat/proteinGroups.txt', quant = "tmt")
rat_pg <- clean_df(rat_pg)

rat_pg_nc <- rat_pg %>% 
  dplyr::select(!(matches("^\\d") & !matches("nc")))

r_pg_long <- rat_pg %>% 
  dplyr::pivot_longer(
    cols = matches("^(reporter|lfq|ratio|intensity)"), 
    names_to = "sample",
    values_to = "log2_intensity"
  ) %>% 
  dplyr::mutate(
    # Remove prefix
    sample = trimws(gsub("(reporter|lfq|ratio|intensity|corrected)", "", sample), which = "left", whitespace = "_"),
    tmt_channel = as.numeric(stringr::str_extract(sample, "^\\d+")),
    sample = gsub("^(\\d+_)", "", sample),
    # Collapse ID columns to first entry
    protein_id_s = gsub("(;.+)", "", protein_id_s),
    gene_names = gsub("(;.+)", "", gene_names),
  ) %>% 
  remove_empty_channels(.) %>% 
  # Log2 transform
  dplyr::mutate(log2_intensity = dplyr::if_else(log2_intensity != 0, log2(log2_intensity), 0))

remove_empty_channels <- function(data){
  ch <- data %>% 
  group_by(tmt_channel) %>% 
  summarise(values = length(unique(log2_intensity))) %>% 
  filter(values <= 2) %>% 
  pull(tmt_channel)
  
  data <- filter(data, !tmt_channel %in% ch)
  
  return(data)
}


write_annotation_file(rat_pg)

protein_groups <- readr::read_tsv("data-raw/proteinGroups.txt",
                                  col_select =  c("Protein IDs",
                                                  "Gene names",
                                                  "Protein names",
                                                  matches("Intensity .*"),
                                                  "Only identified by site",
                                                  "Reverse",
                                                  "Potential contaminant"
                                  ))

proteingroups_clean <- clean_df(protein_groups)

pg_long <- proteingroups_clean %>% 
  select(!only_identified_by_site) %>% 
  pivot_longer(cols = starts_with("intensity"),
               names_to = "experiment",
               values_to = "intensity") %>% 
  mutate(intensity = log2(intensity)) %>% 
  mutate(experiment = str_remove(experiment, "intensity_"))

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
