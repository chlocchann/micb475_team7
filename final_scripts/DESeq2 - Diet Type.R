#### Loading packages #### 
library(DESeq2)
library(phyloseq) 
library(tidyverse)
library(dplyr)
library(tidyr)  # For reshaping data


#### Load phyloseq object ####
load('zoo_rare_10000.RData') 


# Transform counts to add a pseudocount of 1
zoo_plus1 <- transform_sample_counts(zoo_rare, function(x) x + 1)



################### I. DESeq2 CODING + GRAPHING -----------------------------------------------------------------

############# A) HERBIVORES
# Subset samples to herbivores
herbivore_samples <- subset_samples(zoo_plus1, OchHMC == 'H') 

#### DESeq for Captive vs Wild for Herbivores ####
# Create DESeq2 object
deseq_herbivores <- phyloseq_to_deseq2(herbivore_samples, ~ captive_wild)
DESEQ_herbivores <- DESeq(deseq_herbivores)

# Results
res_herbivores <- results(DESEQ_herbivores, tidy=TRUE, 
                          contrast = c("captive_wild", "captive", "wild"))
View(res_herbivores)

# Table of results for significant ASVs
sigASVs_herbivores <- res_herbivores |>
  filter(padj < 0.05 & abs(log2FoldChange) > 1.5) |>
  dplyr::rename(ASV = row)
View(sigASVs_herbivores) ##ensure we do have groups with significant difference btw captive & wild

# Vector of significant ASVs
sigASVs_herbivores_vec <- sigASVs_herbivores %>%
  pull(ASV)

# Identifying the differentially abundant ASVs at the phylum level
herbivores_DESeq <- prune_taxa(sigASVs_herbivores_vec, herbivore_samples)
herbivores_sigASVs <- tax_table(herbivores_DESeq) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_herbivores, by = "ASV") %>%
  arrange(log2FoldChange) %>%
  # Remove 'p__' prefix and any numbers after the phylum name
  mutate(Phylum_Clean = gsub("p__|\\..*", "", Phylum)) %>%
  mutate(Phylum_Clean = factor(Phylum_Clean, levels = unique(Phylum_Clean)))

# Plot at the phylum level
herbivores_DESeq_plot <- ggplot(herbivores_sigASVs) +
  geom_point(aes(x = Phylum_Clean, y = log2FoldChange)) +
  geom_errorbar(aes(x = Phylum_Clean, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Phylum", y = "Log2 Fold Change (captive/wild)") +
  ggtitle("Herbivores (Diet Type)")
herbivores_DESeq_plot

# Export the phylum plot
ggsave(filename = "DESeq_phylum_plot_herbivores.png", herbivores_DESeq_plot,
       height = 5, width = 7)



############# -----------------------------------------------------------

############# B) OMNIVORES
# Subset samples to omnivores
omnivore_samples <- subset_samples(zoo_plus1, OchHMC == 'O')

#### DESeq for Captive vs Wild for Omnivores ####
# Create DESeq2 object
deseq_omnivores <- phyloseq_to_deseq2(omnivore_samples, ~ captive_wild)
DESEQ_omnivores <- DESeq(deseq_omnivores)

# Results
res_omnivores <- results(DESEQ_omnivores, tidy=TRUE, 
                         contrast = c("captive_wild", "captive", "wild"))
View(res_omnivores)

# Table of results for significant ASVs
sigASVs_omnivores <- res_omnivores |>
  filter(padj < 0.05 & abs(log2FoldChange) > 1.5) |>
  dplyr::rename(ASV = row)
View(sigASVs_omnivores) ##ensure we do have groups with significant difference btw captive & wild

# Vector of significant ASVs
sigASVs_omnivores_vec <- sigASVs_omnivores %>%
  pull(ASV)

# Identifying the differentially abundant ASVs at the phylum level
omnivores_DESeq <- prune_taxa(sigASVs_omnivores_vec, omnivore_samples)
omnivores_sigASVs <- tax_table(omnivores_DESeq) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_omnivores, by = "ASV") %>%
  arrange(log2FoldChange) %>%
  # Remove 'p__' prefix and any numbers after the phylum name
  mutate(Phylum_Clean = gsub("p__|\\..*", "", Phylum)) %>%
  mutate(Phylum_Clean = factor(Phylum_Clean, levels = unique(Phylum_Clean)))

# Plot at the phylum level
omnivores_DESeq_plot <- ggplot(omnivores_sigASVs) +
  geom_point(aes(x = Phylum_Clean, y = log2FoldChange)) +
  geom_errorbar(aes(x = Phylum_Clean, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Phylum", y = "Log2 Fold Change (captive/wild)") +
  ggtitle("Omnivores (Diet Type)")
omnivores_DESeq_plot

# Export the phylum plot
ggsave(filename = "DESeq_phylum_plot_omnivores.png", omnivores_DESeq_plot,
       height = 5, width = 7)





############# -----------------------------------------------------------

############# C) CARNIVORES
# Subset samples to carnivores
carnivore_samples <- subset_samples(zoo_plus1, OchHMC == 'C') 

#### DESeq for Captive vs Wild for Carnivores ####
# Create DESeq2 object
deseq_carnivores <- phyloseq_to_deseq2(carnivore_samples, ~ captive_wild)
DESEQ_carnivores <- DESeq(deseq_carnivores)

# Results
res_carnivores <- results(DESEQ_carnivores, tidy=TRUE, 
                          contrast = c("captive_wild", "captive", "wild"))
View(res_carnivores)

# Table of results for significant ASVs
sigASVs_carnivores <- res_carnivores |>
  filter(padj < 0.05 & abs(log2FoldChange) > 1.5) |>
  dplyr::rename(ASV = row)
View(sigASVs_carnivores) #ensure we do have groups with significant difference btw captive & wild

# Vector of significant ASVs
sigASVs_carnivores_vec <- sigASVs_carnivores %>%
  pull(ASV)

# Identifying the differentially abundant ASVs at the phylum level
carnivores_DESeq <- prune_taxa(sigASVs_carnivores_vec, carnivore_samples)
carnivores_sigASVs <- tax_table(carnivores_DESeq) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_carnivores, by = "ASV") %>%
  arrange(log2FoldChange) %>%
  # Remove 'p__' prefix and any numbers after the phylum name
  mutate(Phylum_Clean = gsub("p__|\\..*", "", Phylum)) %>%
  mutate(Phylum_Clean = factor(Phylum_Clean, levels = unique(Phylum_Clean)))

# Plot at the phylum level
carnivores_DESeq_plot <- ggplot(carnivores_sigASVs) +
  geom_point(aes(x = Phylum_Clean, y = log2FoldChange)) +
  geom_errorbar(aes(x = Phylum_Clean, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Phylum", y = "Log2 Fold Change (captive/wild)") +
  ggtitle("Carnivores (Diet Type)")
carnivores_DESeq_plot

# Export the phylum plot
ggsave(filename = "DESeq_phylum_plot_carnivores.png", carnivores_DESeq_plot,
       height = 5, width = 7)







################### II. "POSITIVE EVENTS": noticed some above 0 ASVs for log2FoldChange for categories -------------
              # meaning higher in adaptive than wild
              # which are relatively (unexpected results)



############# A) HERBIVORES
posit_herbivores_ASVs <- herbivores_sigASVs %>%
  filter(log2FoldChange > 0) %>%
  pull(ASV)

posit_herbivores_ASVs
    # [1] "f38b1cb0b64fd29cd983733ecdc2ebe1" "8b10aacaed7b9efb5a000f1c928e2552" "6bcc9245ba0170f675128485fb0bc376"
    # [4] "ca3a12f5eaccca9aa25853a1ab844f63" "aeea90c86c8bf9626cc666181339e2c9" "40420155a410c1d17a505731cd700146"
    # [7] "a98be5b3117ef2fa81141b4d214529a3" "1b00280ce9cd074cd6204adb083fb021" "7cf74cb58f6a831d55bdc8ffaef348ac"
    # [10] "bc511f2c8807bf72a5dd13d5a8dfac1d" "f001e292e9e19ce85465aaf6413bde55" "8e6a4d9e0947a5d5c14fa6c5d5451dd5"
    # [13] "f50defe8cf4a491bc5377919cec7e331" "37afc484bc43bc072162eeb576c474f9" "60f1d2183c3663fc8d385631d75616ff"
    # [16] "bdbefdad98fb1a9d5b63bea068b4db17"

# Prune the phyloseq object to keep only the specified ASVs
target_ASV_herbivores <- c("f38b1cb0b64fd29cd983733ecdc2ebe1", "8b10aacaed7b9efb5a000f1c928e2552",
                           "6bcc9245ba0170f675128485fb0bc376", "ca3a12f5eaccca9aa25853a1ab844f63",
                           "aeea90c86c8bf9626cc666181339e2c9", "40420155a410c1d17a505731cd700146", 
                           "a98be5b3117ef2fa81141b4d214529a3", "1b00280ce9cd074cd6204adb083fb021",
                           "7cf74cb58f6a831d55bdc8ffaef348ac", "bc511f2c8807bf72a5dd13d5a8dfac1d", 
                           "f001e292e9e19ce85465aaf6413bde55", "8e6a4d9e0947a5d5c14fa6c5d5451dd5", 
                           "f50defe8cf4a491bc5377919cec7e331", "37afc484bc43bc072162eeb576c474f9", 
                           "60f1d2183c3663fc8d385631d75616ff", "bdbefdad98fb1a9d5b63bea068b4db17")

target_ASV_phyloseq_herbivores <- prune_taxa(target_ASV_herbivores, herbivore_samples)

# Get the taxonomic information for the specified ASVs
tax_info_herbivores <- tax_table(target_ASV_phyloseq_herbivores) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")
tax_info_herbivores

# Get the OTU table to see the counts for the specified ASVs
otu_counts_herbivores <- otu_table(target_ASV_phyloseq_herbivores) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")
otu_counts_herbivores

# Combine taxonomic information with OTU counts
combined_info_herbivores <- merge(tax_info_herbivores, otu_counts_herbivores, by = "ASV")
combined_info_herbivores

# Identify the columns for sample counts (assuming they start after the Species column)
count_columns_herbivores <- names(combined_info_herbivores)[-(1:8)]  # Adjust based on the number of taxonomy columns

# Reshape combined_info to have sample IDs as a column
combined_info_long_herbivores <- combined_info_herbivores %>%
  pivot_longer(cols = all_of(count_columns_herbivores), # Only include count columns
               names_to = "SampleID",  # New column for sample IDs
               values_to = "Count") # New column for counts

# Get sample data to understand the context of these ASVs
sample_info_herbivores <- sample_data(target_ASV_phyloseq_herbivores) %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID")

# Merge the reshaped combined_info with sample information and view
combined_results_herbivores <- merge(combined_info_long_herbivores, sample_info_herbivores, by = "SampleID", all.x = TRUE)
View(combined_results_herbivores)

# Save the merged results as a CSV file
write.csv(combined_results_herbivores, "positive_ASV_herbivores.csv", row.names = FALSE)





############# -----------------------------------------------------------

############# C) OMNIVORES
posit_omnivores_ASVs <- omnivores_sigASVs %>%
  filter(log2FoldChange > 0) %>%
  pull(ASV)

posit_omnivores_ASVs
        # [1] "055aec5edbd03daea137c9fe44ec23af" "f001e292e9e19ce85465aaf6413bde55" "ca3a12f5eaccca9aa25853a1ab844f63"
        # [4] "963e0574e367692d870c4b02306193be" "40420155a410c1d17a505731cd700146" "844f46a02a277584eae636be55ac3a91"
        # [7] "a6cd3704b9d621add77f9305332721eb"

# Prune the phyloseq object to keep only the specified ASVs
target_ASV_omnivores <- c("055aec5edbd03daea137c9fe44ec23af", "f001e292e9e19ce85465aaf6413bde55",
                          "ca3a12f5eaccca9aa25853a1ab844f63", "963e0574e367692d870c4b02306193be", 
                          "40420155a410c1d17a505731cd700146", "844f46a02a277584eae636be55ac3a91", 
                          "a6cd3704b9d621add77f9305332721eb")

target_ASV_phyloseq_omnivores <- prune_taxa(target_ASV_omnivores, omnivore_samples)

# Get the taxonomic information for the specified ASVs
tax_info_omnivores <- tax_table(target_ASV_phyloseq_omnivores) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")
tax_info_omnivores

# Get the OTU table to see the counts for the specified ASVs
otu_counts_omnivores <- otu_table(target_ASV_phyloseq_omnivores) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")
otu_counts_omnivores

# Combine taxonomic information with OTU counts
combined_info_omnivores <- merge(tax_info_omnivores, otu_counts_omnivores, by = "ASV")
combined_info_omnivores

# Identify the columns for sample counts (assuming they start after the Species column)
count_columns_omnivores <- names(combined_info_omnivores)[-(1:8)]  # Adjust based on the number of taxonomy columns

# Reshape combined_info to have sample IDs as a column
combined_info_long_omnivores <- combined_info_omnivores %>%
  pivot_longer(cols = all_of(count_columns_omnivores), # Only include count columns
               names_to = "SampleID",  # New column for sample IDs
               values_to = "Count") # New column for counts

# Get sample data to understand the context of these ASVs
sample_info_omnivores <- sample_data(target_ASV_phyloseq_omnivores) %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID")

# Merge the reshaped combined_info with sample information and view
combined_results_omnivores <- merge(combined_info_long_omnivores, sample_info_omnivores, by = "SampleID", all.x = TRUE)
View(combined_results_omnivores)

# Save the merged results as a CSV file
write.csv(combined_results_omnivores, "positive_ASV_omnivores.csv", row.names = FALSE)




############# -----------------------------------------------------------

############# C) CARNIVORES
posit_carnivores_ASVs <- carnivores_sigASVs %>%
  filter(log2FoldChange > 0) %>%
  pull(ASV)

posit_carnivores_ASVs
      # [1] "cc86fc399e282b8bfda8ca18779fffac" "55517091a9f60137000fbb6d18e9dc86" "47085e3c4d82d94df39fe6519aae94d5"
      # [4] "571c5868ce894564acbf7c1e07fd65ff" "d8b9abc78138c51d05f80472e8914f14" "40420155a410c1d17a505731cd700146"
      # [7] "3443833956efe96c2c69cc0d0f761f8f" "bf85fc4661d163e33ec8d844c01515ec" "b545600dffab681ebe77a5024816a9b5"
      # [10] "251d3205b39d1916cb63897d12089b83" "d5366f4e78d8e11e5e3ca2134d48f3bd" "5afef9c3934a3fd690f59a2675b17ad2"
      # [13] "53d42797b8ac1fa460a02c8a2904581e" "f50defe8cf4a491bc5377919cec7e331"

# Prune the phyloseq object to keep only the specified ASVs
target_ASV_carnivores <- c("cc86fc399e282b8bfda8ca18779fffac", "55517091a9f60137000fbb6d18e9dc86", 
                           "47085e3c4d82d94df39fe6519aae94d5", "571c5868ce894564acbf7c1e07fd65ff", 
                           "d8b9abc78138c51d05f80472e8914f14", "40420155a410c1d17a505731cd700146", 
                           "3443833956efe96c2c69cc0d0f761f8f", "bf85fc4661d163e33ec8d844c01515ec", 
                           "b545600dffab681ebe77a5024816a9b5", "251d3205b39d1916cb63897d12089b83", 
                           "d5366f4e78d8e11e5e3ca2134d48f3bd", "5afef9c3934a3fd690f59a2675b17ad2", 
                           "53d42797b8ac1fa460a02c8a2904581e", "f50defe8cf4a491bc5377919cec7e331")

target_ASV_phyloseq_carnivores <- prune_taxa(target_ASV_carnivores, carnivore_samples)

# Get the taxonomic information for the specified ASVs
tax_info_carnivores <- tax_table(target_ASV_phyloseq_carnivores) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")
tax_info_carnivores

# Get the OTU table to see the counts for the specified ASVs
otu_counts_carnivores <- otu_table(target_ASV_phyloseq_carnivores) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")
otu_counts_carnivores

# Combine taxonomic information with OTU counts
combined_info_carnivores <- merge(tax_info_carnivores, otu_counts_carnivores, by = "ASV")
combined_info_carnivores

# Identify the columns for sample counts (assuming they start after the Species column)
count_columns_carnivores <- names(combined_info_carnivores)[-(1:8)]  # Adjust based on the number of taxonomy columns

# Reshape combined_info to have sample IDs as a column
combined_info_long_carnivores <- combined_info_carnivores %>%
  pivot_longer(cols = all_of(count_columns_carnivores), # Only include count columns
               names_to = "SampleID",  # New column for sample IDs
               values_to = "Count") # New column for counts

# Get sample data to understand the context of these ASVs
sample_info_carnivores <- sample_data(target_ASV_phyloseq) %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID")

# Merge the reshaped combined_info with sample information and view
combined_results_carnivores <- merge(combined_info_long_carnivores, sample_info_carnivores, by = "SampleID", all.x = TRUE)
View(combined_results_carnivores)

# Save the merged results as a CSV file
write.csv(combined_results_carnivores, "positive_ASV_carnivores.csv", row.names = FALSE)

