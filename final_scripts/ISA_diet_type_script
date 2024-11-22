library(tidyverse)
library(phyloseq)
library(indicspecies)

load("zoo_rare_10000.RData")

# Aggregate taxa to the species level & transforming to relative abundance
species <- tax_glom(zoo_rare, "Species", NArm = FALSE)
species_RA <- transform_sample_counts(species, function(x) x / sum(x))

# clustering Diet type (OchHMC) and Captivity
metadata_sample_diet <- sample_data(species_RA) %>% as.data.frame()   #extracting the metadata
metadata_sample_diet$combined_cluster <- paste(metadata_sample$OchHMC, metadata_sample$captive_wild, sep = "_")
sample_data(species_RA) <- sample_data(metadata_sample_diet)

# Run the Indicator Species Analysis with the clustered variable (GutFermentersHMC and captivity)
diet_isa_output <- multipatt(t(otu_table(species_RA)), cluster = metadata_sample_diet$combined_cluster)

# Merging taxonomy table with the phyloseq object & filtering by significant p-value
taxtable <- tax_table(zoo_rare) %>% as.data.frame() %>% rownames_to_column(var="ASV")
diet_isa_results <- diet_isa_output$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05, stat >= 0.80) 
write.csv(diet_isa_results, "ISA_DIET_results_unfiltered.csv")

### CAPTIVE CARNIVORE ###
captive_C <- diet_isa_results %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, s.C_captive, p.value, stat) %>%
  filter(s.C_captive != 0) %>% 
  mutate(stat = signif(stat, digits = 3)) %>%
  arrange(Phylum, p.value, stat) %>%
  select(-s.C_captive)

write.csv(captive_C, "ISA_captive_carnivore.csv", row.names = FALSE)

### CAPTIVE HERBIVORE ###
captive_H <- diet_isa_results %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, s.H_captive, p.value, stat) %>%
  filter(s.H_captive != 0) %>% 
  mutate(stat = signif(stat, digits = 3)) %>%
  arrange(Phylum, p.value, stat) %>%
  select(-s.H_captive)

write.csv(captive_H, "ISA_captive_herbivore.csv", row.names = FALSE)

### combining both captive tables ###
combined_captive <- bind_rows(captive_C, captive_H) %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, s.C_captive, s.H_captive, p.value, stat)
write.csv(combined_captive, "ISA_combined_captive.csv", row.names = FALSE)

### WILD CARNIVORE ###
wild_C <- diet_isa_results %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, s.C_wild, p.value, stat) %>%
  filter(s.C_wild != 0) %>%
  mutate(stat = signif(stat, digits = 3)) %>%
  arrange(Phylum, p.value, stat)  %>%
  select(-s.C_wild)

write.csv(wild_C, "ISA_wild_C.csv", row.names = FALSE)

### WILD HERBIVORE ###
wild_H <- diet_isa_results %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, s.H_wild, p.value, stat) %>%
  filter(s.H_wild != 0) %>%
  mutate(stat = signif(stat, digits = 3)) %>%
  arrange(Phylum, p.value, stat) %>%
  select(-s.H_wild)

write.csv(wild_H, "ISA_wild_H.csv", row.names = FALSE)

### combining both wild tables ###
combined_wild <- bind_rows(wild_C, wild_H) %>%
  select(ASV, Phylum, Class, Order, Family, Genus, Species, s.C_wild, s.H_wild, p.value, stat)
write.csv(combined_wild, "ISA_combined_wild.csv", row.names = FALSE)
