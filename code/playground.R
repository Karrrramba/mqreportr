library(ggthemes)
library(ggokabeito)
library(ggdist)
library(patchwork)


rat_pg_annotated <- load_annotations(a_file = "annotation_file.csv", rat_pg_long) 

# proteinGroups----
rat_pg <- load_table('data-raw/tmt_rat/proteinGroups.txt', quant = "tmt")
rat_pg <- clean_df(rat_pg)
rat_pg_long <- transform_table(rat_pg)

rat_pg_nc <- rat_pg_long %>% 
  filter(grepl("nc", sample, fixed = TRUE))

write_annotation_file(rat_pg_nc)

rat_pg_annotated <- load_annotations(annotation = "annotation_file.csv", data = rat_pg_nc)


# NEM sites----
rat_nem <- load_table("data-raw/tmt_rat/NEM (C)Sites.txt", quant = 'tmt')
rat_nem <- clean_df(rat_nem, prob_treshold = 0.75)
rat_nem_long <- transform_long(rat_nem, site = TRUE)

rat_nem_long <- rat_nem_long %>% 
  filter(!sample %in% as.character(seq(1:11)))

write_group_annotations_file(rat_nem_long)
  
rat_nem_groups <- load_group_annotations("rat_nem_long_annotations_2024-08-09 20.09.csv", rat_nem_long)

rat_nem_batched <- rat_nem_groups %>%
  lazy_dt() %>% 
  mutate(batch = str_extract(sample, "([a-z0-9]+)$")) %>% 
  filter(full_name != "empty") %>% 
  as.data.table()
sorted <- rat_nem_batched %>% 
  lazy_dt() %>% 
  group_by(batch) %>% 
  arrange(batch) %>%
  ungroup() %>% 
  as.data.table()

# box plots
# non- normalized
box_none <- rat_nem_batched %>% 
  as.data.table() %>% 
  ggplot(aes(x = interaction(full_name, batch), y = log2_intensity, color = batch)) +
  geom_boxplot() +
  theme_tufte() +
  ggtitle("Median-normalized data")

# median-normalized
box_median <- rat_nem_batched %>% 
  lazy_dt() %>% 
  group_by(full_name) %>% 
  mutate(log2_intensity = log2_intensity -  median(log2_intensity, na.rm = TRUE)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  ggplot(aes(x = interaction(full_name, batch), y = log2_intensity, color = batch)) +
  geom_boxplot() +
  theme_tufte() +
  ggtitle("Median-normalized data")

box_none + box_median + plot_layout(nrow = 2, axes = "collect", axis_titles = "collect") + 

rat_nem_batched %>% 
  lazy_dt() %>% 
  group_by(full_name) %>% 
  arrange(batch) %>% 
  mutate(log2_intensity = log2_intensity -  median(log2_intensity, na.rm = TRUE)) %>% 
  ungroup() %>% 
  as.data.table() %>% 
  ggplot(aes(x = interaction(full_name, batch), y = log2_intensity, color = batch)) +
  geom_boxplot() +
  theme_tufte() +
  ggtitle("Median-normalized data")


rat_nem_wide <- rat_nem_groups %>% 
  lazy_dt() %>% 
  filter(!is.na(age) & !str_detect(full_name, "CF")) %>% #remove CF
  mutate(
    age = gsub("\\D+", "", age),
    sample_name = paste0(sex, "_", diet, "_", age)
    ) %>% 
  select(!c(sample, sex, diet, age, channel)) %>% 
  rename("sample" = "full_name") %>% 
  pivot_wider(
    id_cols = 1:7,
    names_from = sample,
    values_from = log2_intensity
  ) %>% 
  as.data.table()



  
  
  
