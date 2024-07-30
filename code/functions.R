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

load_table <- function(df, multi) {
  
  cols <- c("Protein IDs",
            "Reverse",
            "Potential contaminant")
  if (stringr::str_detect(df, "site")) {
    cols <- append(cols, c())
  }
  
  if (multi == "tmt") {
    quant_columns <- paste0("matches(", "^(LFQ) |(I|i)ntensity .*", ")")
    cols <- append(cols, quant_columns)
  } else if (multi == "lfq") {
    quant_columns <- paste0("matches(", "Reporter intensity (corrected|\\d)", ")")
    cols <- append(cols, quant_columns)
  } else {
    quant_columns <- paste0("matches(", "Ratio .*", ")")
    cols <- append(cols, quant_columns)
  }
  
  cols <- cols[!grepl("count", cols)]
  
  df <- readr::read_tsv(file = df,
                        col_select = cols)
  return(df)
}

list2env(
  lapply(setNames(site_tables, janitor::make_clean_names(site_tables)), load_table, multi = "lfq"),
  envir = .GlobalEnv)

### Universal data cleaning
# Formats column names, removes reverse and contaminants 
clean_df <- function(data) {
  del_row <- which(data[, "Reverse"] == "+" | data[, "Potential contaminant"] == "+")
  del_col <- which(names(data) %in% c("Reverse", "Potential contaminant"))
  
  data <- data[-del_row, -del_col]
  data <- janitor::clean_names(data, abbreviation = "ID")
  return(data)
}
# proteinGroups
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
