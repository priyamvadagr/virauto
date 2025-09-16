##########################################################################################################################
#   VirScan Data Analysis – Full peptide-level output (with metadata)
#   Adapted from Cortese et al. (2021) EBV-MS manuscript
#
#   Modifications:
#   - Includes ALL peptides (no filtering by virus or significance)
#   - Each peptide (pep_id) = one row
#   - Metadata columns joined for interpretation
##########################################################################################################################

library(tidyverse)
library(broom)
library(writexl)

######################################################
# 1. Load dataset and create metadata
######################################################
setwd('/ix/djishnu/Priyamvada/virauto/data/epitopes')
VirScan_dataset <- read_csv("Virscan_dataset.csv")

Metadata <- VirScan_dataset %>% 
  select(-id1, -casestat, -timeser, -zscore, -repnum) %>% 
  group_by(pep_id) %>% 
  slice_head(n = 1) %>% 
  mutate(pep_id = as.character(pep_id))

#################################################################################################################
# 2. Fisher’s Exact Test: case-control difference in positive hits
#################################################################################################################

########################### PRE-ONSET SAMPLES ###########################

epihit_bothreps_preonset <- VirScan_dataset %>% 
  filter(timeser == "last") %>%  
  select(casestat, id1, pep_id, repnum, timeser, zscore) %>% 
  arrange(pep_id, id1) %>% 
  mutate(
    epi_hit = if_else(zscore > 3.5, 1, 0)
  ) %>% 
  group_by(id1, pep_id) %>% 
  mutate(hit_both = if_else(sum(epi_hit) == 2, 1, 0))

hit_in_both_pre <- epihit_bothreps_preonset %>% 
  ungroup() %>% 
  filter(repnum == 1) %>%  # keep only one replicate row
  group_by(pep_id) %>% 
  filter(sum(hit_both) >= 1)

hit_in_both_final_pre <- hit_in_both_pre %>% 
  group_by(pep_id, casestat) %>% 
  summarise(prop_with_hit = mean(hit_both >= 1), .groups = "drop") %>% 
  pivot_wider(names_from = casestat, values_from = prop_with_hit) %>% 
  rename(Prop_controls = "0", Prop_cases = "1") %>% 
  mutate(
    Cases_with_hit = Prop_cases * 30,
    Cases_without_hit = 30 - Cases_with_hit,
    Controls_with_hit = Prop_controls * 30,
    Controls_without_hit = 30 - Controls_with_hit
  )

# Fisher’s exact test
hit_in_both_final_fisher_pre <- hit_in_both_final_pre %>% 
  group_by(pep_id) %>% 
  nest() %>% 
  mutate(matrix = map(data, ~ matrix(unlist(.x), nrow = 2))) %>% 
  mutate(fisher = map(matrix, ~ fisher.test(.x))) %>% 
  mutate(stats = map(fisher, ~ broom::tidy(.x))) %>% 
  unnest(c(data, stats)) %>% 
  select(pep_id, Cases_with_hit, Controls_with_hit, p.value) %>% 
  mutate(
    Prop_cases = Cases_with_hit / 30,
    Prop_controls = Controls_with_hit / 30,
    High_in_cases = if_else(Prop_cases > Prop_controls, 1, 0),
    pep_id = as.character(pep_id)
  ) %>% 
  ungroup()

# Add metadata
full_pre_table <- hit_in_both_final_fisher_pre %>% 
  left_join(Metadata, by = "pep_id")

write_xlsx(full_pre_table, "pre_onset_case_control_compare.xlsx")

########################### POST-ONSET SAMPLES ###########################

epihit_bothreps_postonset <- VirScan_dataset %>% 
  filter(timeser == "after") %>%  
  select(casestat, id1, pep_id, repnum, timeser, zscore) %>% 
  arrange(pep_id, id1) %>% 
  mutate(
    epi_hit = if_else(zscore > 3.5, 1, 0)
  ) %>% 
  group_by(id1, pep_id) %>% 
  mutate(hit_both = if_else(sum(epi_hit) == 2, 1, 0))

hit_in_both_post <- epihit_bothreps_postonset %>% 
  ungroup() %>% 
  filter(repnum == 1) %>% 
  group_by(pep_id) %>% 
  filter(sum(hit_both) >= 1)

hit_in_both_final_post <- hit_in_both_post %>% 
  group_by(pep_id, casestat) %>% 
  summarise(prop_with_hit = mean(hit_both >= 1), .groups = "drop") %>% 
  pivot_wider(names_from = casestat, values_from = prop_with_hit) %>% 
  rename(Prop_controls = "0", Prop_cases = "1") %>% 
  mutate(
    Cases_with_hit = Prop_cases * 30,
    Cases_without_hit = 30 - Cases_with_hit,
    Controls_with_hit = Prop_controls * 30,
    Controls_without_hit = 30 - Controls_with_hit
  )

# Fisher’s exact test
hit_in_both_final_fisher_post <- hit_in_both_final_post %>% 
  group_by(pep_id) %>% 
  nest() %>% 
  mutate(matrix = map(data, ~ matrix(unlist(.x), nrow = 2))) %>% 
  mutate(fisher = map(matrix, ~ fisher.test(.x))) %>% 
  mutate(stats = map(fisher, ~ broom::tidy(.x))) %>% 
  unnest(c(data, stats)) %>% 
  select(pep_id, Cases_with_hit, Controls_with_hit, p.value) %>% 
  mutate(
    Prop_cases = Cases_with_hit / 30,
    Prop_controls = Controls_with_hit / 30,
    High_in_cases = if_else(Prop_cases > Prop_controls, 1, 0),
    pep_id = as.character(pep_id)
  ) %>% 
  ungroup()

# Add metadata
full_post_table <- hit_in_both_final_fisher_post %>% 
  left_join(Metadata, by = "pep_id")

write_xlsx(full_post_table, "post_onset_case_control_compare.xlsx")
