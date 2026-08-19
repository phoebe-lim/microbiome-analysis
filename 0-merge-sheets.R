library(tidyverse)
library(janitor)
setwd("~/OneDrive - Imperial College London/MSc Genomic Medicine/Research Project/Microbiome Analysis")

# Load files --------------------------------------------------------------
pf <- read_csv("PROFOUND/metadata_PFND.csv")
pf_map <- read_tsv("PROFOUND/Manifest-Map_BAL.txt")
bru <- read_csv("BRU/metadata_BRU.csv")
master_sheet <- read_csv("Metadata Processing/master_sheet.csv")

# Merge PROFOUND metadata and manifest map --------------------------------
# Only keep relevant columns in manifest map
pf_map <- pf_map %>%
  select(c("sample-id", "DNA (ng/ul)", "SampleType", "Diagnosis")) %>%
  rename("sample_id" = "sample-id",
         "ddPCR" = "DNA (ng/ul)",
         "sample_type" = "SampleType",
         "diagnosis" = "Diagnosis")

pf <- pf %>%
  mutate(patient_id = str_c(patient_id, "-BAL")) %>%
  select(-diagnosis)

# Keep all samples in manifest map
pf_merged <- pf_map %>%
  left_join(pf, by = c("sample_id" = "patient_id"))

# Merge PROFOUND and BRU metadata -----------------------------------------
colnames(pf_merged)
colnames(bru)

# Align PROFOUND columns with BRU
pf_merged <- pf_merged %>%
  mutate(death_1y = if_else(death == "Yes" & survival_days < 366,
                                "Yes", "No")) %>%
  rename(MUC5B = genotype) %>%
  select(-c(FVC_1y, ppFVC_1y, DLCO_1y, ppDLCO_1y, death)) # drop columns not in BRU

# Drop BRU columns not in PROFOUND and match column names
bru <- bru %>%
  select(-c(Baseline_FEV1, Baseline_ppFEV1)) %>%
  rename("sample_id" = "sample-id",
         diagnosis = Diagnosis,
         sample_type = Description,
         age = Age,
         sex = Sex,
         smoking_hist = Smoking_history,
         FVC_baseline = Baseline_FVC,
         ppFVC_baseline = Baseline_ppFVC,
         DLCO_baseline = Baseline_TLCO,
         ppDLCO_baseline = Baseline_ppTLCO,
         death_1y = Death_1y,
         status_1y = Status_1y,
         survival_days = Survival_days)

# Standardise values in sample_type and diagnosis columns
bru <- bru %>%
  mutate(sample_type = str_replace_all(sample_type, "^BKG.*", "BKG") %>% # remove all text after "BKG"
           str_replace_all("\\.", " ") %>% # replace all full stops with spaces
           replace_values("BAL Flush Control" ~ "BAL Flush",
                          "Water oral control" ~ "Water Oral Control"),
         diagnosis = replace_values(diagnosis, "Negative" ~ "Control"))

# Add column indicating run
pf_merged <- mutate(pf_merged, run = "PFND")
bru <- mutate(bru, run = "BRU")

pf_bru <- merge(pf_merged, bru, all = TRUE) %>%
  distinct() # drop duplicate rows
write_csv(pf_bru, "Metadata Processing/metadata_BRU_PFND.csv")

# Format combined BRU, PROFOUND and PROFILE metadata sheet --------------------
colnames(master_sheet)

# Change column names
master_sheet <- master_sheet %>%
  select(-c(Baseline_FEV1, Baseline_ppFEV1)) %>%
  rename("sample_id" = "sample-id",
         diagnosis = Diagnosis,
         sample_type = Description,
         age = Age,
         sex = Sex,
         smoking_history = Smoking_history,
         FVC_baseline = Baseline_FVC,
         ppFVC_baseline = Baseline_ppFVC,
         DLCO_baseline = Baseline_TLCO,
         ppDLCO_baseline = Baseline_ppTLCO,
         death_1y = Death_1y,
         status_1y = Status_1y,
         survival_days = Survival_days)

# Check unique values of categorical columns and their counts
master_sheet %>%
  select(diagnosis, sample_type, smoking_history) %>% 
  map(table)

# Correct typos and standardise values
master_sheet <- master_sheet %>%
  mutate(diagnosis = replace_when(diagnosis,
                                  str_detect(diagnosis, regex("IPF", ignore_case = TRUE)) ~ "IPF") %>% 
                                  # drop additional diagnoses on top of IPF
           replace_values("Familai fNSIP" ~ "NSIP"),
         sample_type = str_replace_all(sample_type, "^BKG.*", "BKG") %>% # remove all text after "BKG"
           str_replace_all("\\.", " ") %>% # replace all full stops with spaces
           replace_values("IPF" ~ "BAL",
                          "Water oral control" ~ "Water Oral Control"),
         smoking_history = str_remove_all(smoking_history, ".smoker"))

# Check for duplicates
get_dupes(master_sheet)
master_sheet <- distinct(master_sheet)
write_csv(master_sheet, "Metadata Processing/metadata_final.csv")

# Fill in missing ddPCR values using BRU + PROFOUND metadata sheet
missing <- read_csv("Metadata Processing/metadata_final.csv")
pf_bru <- read_csv("Metadata Processing/metadata_BRU_PFND.csv")

filled <- missing %>%
  left_join(select(pf_bru, c(sample_id, ddPCR)), by = "sample_id", suffix = c("", "_pf_bru")) %>%
  mutate(ddPCR = coalesce(ddPCR, ddPCR_pf_bru)) %>%
  select(-ddPCR_pf_bru)

# Add column indicating study
filled <- mutate(filled, study = case_when(sample_type == "BAL" & str_starts(key, "BRU")  ~ "BRU",
                                           sample_type == "BAL" & str_starts(key, "PFND") ~ "PROFOUND",
                                           sample_type == "BAL" & str_starts(key, "PRO")  ~ "PROFILE",
                                           .default = "Control"))

write_csv(filled, "Metadata Processing/metadata_final.csv")
