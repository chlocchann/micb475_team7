#!/usr/bin/env Rscript

#### Installing and loading packages #### 
library(tidyverse)
library(dplyr)
library(readr)

#### Load Data ####
zoo_meta <- read_delim("zoo_metadata.tsv")

#### Wrangling Metadata ####
# rename `21-1_PopulationDensity_n/km2` to `21-1_PopulationDensity_n_km2`
colnames(zoo_meta)[55]  <- "21-1_PopulationDensity_n_km2"

# rename `27-1_HuPopDen_Min_n/km2` to `27-1_HuPopDen_Min_n_km2`
colnames(zoo_meta)[73]  <- "27-1_HuPopDen_Min_n_km2"

# rename `27-2_HuPopDen_Mean_n/km2` to `27-3_HuPopDen_5p_n/km2`
colnames(zoo_meta)[74]  <- "27-2_HuPopDen_Mean_n_km2"

# rename `27-3_HuPopDen_5p_n/km2` to `27-3_HuPopDen_5p_n_km2`
colnames(zoo_meta)[75]  <- "27-3_HuPopDen_5p_n_km2"

# append captivity status to OchHMC column
zoo_meta <- zoo_meta %>%
  mutate(
    OchHMC = case_when(
      captive_wild == "captive" ~ paste(OchHMC, "- CAPTIVE"),
      captive_wild == "wild" ~ paste(OchHMC, "- WILD")
    )
  )

# append captivity status to GutFermentersHMC column
zoo_meta <- zoo_meta %>%
  mutate(
    GutFermentersHMC = case_when(
      captive_wild == "captive" ~ paste(GutFermentersHMC, "- CAPTIVE"),
      captive_wild == "wild" ~ paste(GutFermentersHMC, "- WILD")
    )
  )

#### Downloading as a .csv ####
write.csv(zoo_meta,"~/Downloads/zoo_metadata.csv", row.names = FALSE)
