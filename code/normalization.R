library(COMBAT) 
library(limma) #loess
library(mixOmics) #NIPALS PCA
library(nipals) #NIPALS PCA
library(harmonizer)
library(missForest)
library(tidyverse)
library(readxl)

benchmark_data <- read_excel("data-raw/benchmark_datasets/PXD013277/TMT_EColi_spike_in_to_ MCF7_proteingroups.xlsx", .name_repair = janitor::make_clean_names)

bench_data_clean <- benchmark_data %>% 
  select(protein_accession, starts_with("TMT"), -ends_with("ms")) %>% 
  mutate(organism = str_extract(protein_accession, "([A-Z]+)$")) %>%
  pivot_longer(cols = starts_with("tmt"),
               names_to = "sample",
               values_to = "intensity") %>% 
  mutate(intensity = log2(intensity),
         spike_in = case_when(
           str_ends(sample, "126") | str_ends(sample, "127n") | str_ends(sample, "127c") ~ 7,
           str_ends(sample, "128n") | str_ends(sample, "128c") | str_ends(sample, "129n") | str_ends(sample, "129c") ~ 15,
           str_ends(sample, "130n") | str_ends(sample, "130c") | str_ends(sample, "131") ~ 45
         ))

# boxplots
ggplot(bench_data_clean, aes(x = interaction(sample, spike_in), y = intensity, colour = spike_in, fill = organism)) +
  geom_boxplot() +
  ylab("Intensity") +
  xlab("Sample") +
  coord_flip()

# Initial NIPALS PCA
# wide and matrix
data_wide <- bench_data_clean %>% 
  select(protein_accession, sample, intensity) %>% 
  pivot_wider(names_from = "sample",
              values_from = "intensity") %>% 
  mutate(across(starts_with("tmt"), ~as.numeric(.))) 
matrix_init <- data.matrix(data_wide)
# set protein accession as row names
row.names(matrix_init) <- matrix_init[, 1]
matrix_init <- data_matrix[, -1]

pca_init <- nipals::nipals(matrix_init, fitted = TRUE)
biplot(pca_init$scores, pca_init$loadings)

# median normalization
data_mean <- data_wide %>% 
  mutate(across(starts_with("tmt"), ~ .x - mean(.x, na.rm = TRUE)))
matrix_mean <- data.matrix(data_mean)

  
data_median <- data_wide %>% 
  mutate(across(starts_with("tmt"), ~ .x - median(.x, na.rm = TRUE)))
matrix_median <- data.matrix(data_median)


data_harm

# amputation