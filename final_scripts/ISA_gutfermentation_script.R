library(tidyverse)
library(phyloseq)
library(indicspecies)

load("zoo_rare_10000.RData")

# Aggregate taxa to the species level & transforming to relative abundance
species <- tax_glom(zoo_rare, "Species", NArm = FALSE)
species_RA <- transform_sample_counts(species, function(x) x / sum(x))

########## GUT FERMENTATION ISA #####################

# clustering GutFermentersHMC and Captivity
metadata_sample <- sample_data(species_RA) %>% as.data.frame()   #extracting the metadata
metadata_sample$combined_cluster <- paste(metadata_sample$GutFermentersHMC, metadata_sample$captive_wild, sep = "_")
sample_data(species_RA) <- sample_data(metadata_sample)

# Run the Indicator Species Analysis with the clustered variable (GutFermentersHMC and captivity)
gut_isa_output <- multipatt(t(otu_table(species_RA)), cluster = metadata_sample$combined_cluster)

# Merging taxonomy table with the phyloseq object & filtering by significant p-value
taxtable <- tax_table(zoo_rare) %>% as.data.frame() %>% rownames_to_column(var="ASV")
gut_isa_results <- gut_isa_output$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05, stat >= 0.80) 
write.csv(isa_results, "ISA_GUT_results_unfiltered.csv")

### HINDGUT CAPTIVE ###
hindgut_captive <- gut_isa_results %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, s.Hindgut_captive, p.value, stat) %>%
  filter(s.Hindgut_captive != 0) %>% 
  mutate(stat = signif(stat, digits = 3)) %>%
  arrange(Phylum, p.value, stat) %>%
  select(-s.Hindgut_captive)
  
write.csv(hindgut_captive, "ISA_hindgut_captive.csv", row.names = FALSE)

### HINDGUT WILD ###
hindgut_wild <- gut_isa_results %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, s.Hindgut_wild, p.value, stat) %>%
  filter(s.Hindgut_wild != 0) %>% 
  mutate(stat = signif(stat, digits = 3)) %>%
  arrange(Phylum, p.value, stat) %>%
  select(-s.Hindgut_wild)

write.csv(hindgut_wild, "ISA_hindgut_wild.csv", row.names = FALSE)

### combining both hindgut wild and captive tables ###
combined_hindgut <- bind_rows(hindgut_captive, hindgut_wild) %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, s.Hindgut_captive, s.Hindgut_wild, p.value, stat)
write.csv(combined_hindgut, "ISA_combined_hindgut.csv", row.names = FALSE)

### FOREGUT CAPTIVE ###
foregut_captive <- gut_isa_results %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, s.Foregut_captive, p.value, stat) %>%
  filter(s.Foregut_captive != 0) %>%
  mutate(stat = signif(stat, digits = 3)) %>%
  arrange(Phylum, p.value, stat) %>%
  select(-s.Foregut_captive)

write.csv(foregut_captive, "ISA_foregut_captive.csv", row.names = FALSE)

### FOREGUT WILD ###
foregut_wild <- gut_isa_results %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, s.Foregut_wild, p.value, stat) %>%
  filter(s.Foregut_wild != 0) %>%
  mutate(stat = signif(stat, digits = 3)) %>%
  arrange(Phylum, p.value, stat) %>%
  select(-s.Foregut_wild)

write.csv(foregut_wild, "ISA_foregut_wild.csv", row.names = FALSE)

### combining both foregut wild and captive tables ###
combined_foregut <- bind_rows(foregut_captive, foregut_wild) %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, s.Foregut_captive, s.Foregut_wild, p.value, stat)
write.csv(combined_foregut, "ISA_combined_foregut.csv", row.names = FALSE)

### NON-FERMENTERS CAPTIVE ###
nonfermenters_captive <- gut_isa_results %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, s.N_captive, p.value, stat) %>%
  filter(s.N_captive != 0) %>% 
  mutate(stat = signif(stat, digits = 3)) %>%
  arrange(Phylum, p.value, stat) %>%
  select(-s.N_captive)

write.csv(nonfermenters_captive, "ISA_nonfermenters_captive.csv", row.names = FALSE)

### NON-FERMENTERS WILD ###
nonfermenters_wild <- gut_isa_results %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, s.N_wild, p.value, stat) %>%
  filter(s.N_wild != 0) %>% 
  mutate(stat = signif(stat, digits = 3)) %>%
  arrange(Phylum, p.value, stat) %>%
  select(-s.N_wild)

write.csv(nonfermenters_wild, "ISA_nonfermenters_wild.csv", row.names = FALSE)

### combining both nonfermenters wild and captive tables ###
combined_nonfermenter <- bind_rows(nonfermenters_captive, nonfermenters_wild) %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, s.N_captive, s.N_wild, p.value, stat)
write.csv(combined_nonfermenter, "ISA_combined_nonfermenter.csv", row.names = FALSE)
