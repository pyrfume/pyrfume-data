if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyr, dplyr, readxl, webchem)

Keller_2016 <- read_excel("12868_2016_287_MOESM1_ESM.xlsx", skip = 2)

subjects <- Keller_2016 %>%
  rename(Subjects = "Subject # (this study)") %>%
  rename(Race = 9) %>%
  mutate(Experiment = "keller2016") %>%
  select(Subjects, Experiment, Gender, Race, Ethnicity, Age) %>%
  distinct(Subjects, .keep_all = TRUE)

write.csv(subjects, "subjects.csv")


stimuli <- Keller_2016 %>% 
  rename(Ratio = 5) %>%
  mutate(Ratio = gsub(",", "", Ratio)) %>%
  select(CIDs = CID, Ratio) %>%
  distinct(CIDs, Ratio) %>%
  mutate(Concentration = sapply(Ratio, function(x) {
    # Split the fraction into numerator and denominator
    parts <- strsplit(x, "/")[[1]]
    # Convert to numeric and divide
    as.numeric(parts[1]) / as.numeric(parts[2])
  })) %>%
  mutate(Solvent = "paraffin oil") %>%
  mutate(Experiment = "keller2016") %>%
  mutate(Stimulus = row_number()) %>%
  select(Stimulus, Experiment, CIDs, Concentration, Ratio, Solvent)
  
write.csv(stimuli, "stimuli.csv")


molecules <- Keller_2016 %>%
  select(CID=3, OdorName=4, CAS=1) %>%
  distinct() %>%
  mutate(CID = as.numeric(CID)) %>%
  mutate_at(vars(CID), as.character) %>%
  left_join(
    pc_prop(.$CID)[c("CID", "CanonicalSMILES", "MolecularWeight")] %>% 
      mutate_at(vars(CID), as.character),
    by = "CID"
  ) %>%
  mutate(OdorName = gsub(" \\(replicate\\)", "", OdorName)) %>%
  distinct() %>%
  mutate(CID = ifelse(is.na(CID) & CAS == "109-19-0", 8038, CID)) %>%
  mutate(CID = ifelse(is.na(CID) & CAS == "3796-70-1", 1549778, CID))

write.csv(molecules, "molecules.csv")


behavior <- Keller_2016 %>%
  rename(Ratio = 5) %>%
  mutate(Ratio = gsub(",", "", Ratio)) %>%
  rename(Subject = "Subject # (this study)") %>%
  rename(CIDs = CID) %>%
  mutate(Experiment = "keller2016") %>%
  select(-1, -2, -4, -7, -8, -9, -10, -11, -12) %>%
  left_join(stimuli, by = c("CIDs", "Ratio")) %>%
  select(-1, -2, -Experiment.y, -Concentration, -Solvent) %>%
  select(Stimulus, Subject, Experiment = Experiment.x, everything()) %>%
  mutate(across(-c(Stimulus, Subject, Experiment), as.character)) %>%
  pivot_longer(
    cols = -c(Stimulus, Subject, Experiment),  # Exclude the first 3 columns
    names_to = "MeasurementValue",  # Name of the new key column
    values_to = "Value"  # Name of the new value column
  )
  
write.csv(behavior, "behavior.csv")

