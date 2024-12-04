#### Loading packages #### 
library(DESeq2)
library(phyloseq) 
library(tidyverse)
library(dplyr)
library(tidyr)  # For reshaping data
library(ggplot2)

#### Load phyloseq object ####
load('zoo_rare_10000.RData') 


# Transform counts to add a pseudocount of 1
zoo_plus1 <- transform_sample_counts(zoo_rare, function(x) x + 1)



#--------------------------- I. DESeq2 CODING + GRAPHING 

############# A) HERBIVORES
# Subset samples to herbivores
herbivores_samples <- subset_samples(zoo_plus1, OchHMC == 'H')

# Subset samples to herbivores
herbivores_plus1 <- transform_sample_counts(herbivores_samples, function(x) x + 1)

# Create DESeq2 object
herbivores_deseq <- phyloseq_to_deseq2(herbivores_plus1, ~ captive_wild)
DESEQ_herbivores <- DESeq(herbivores_deseq)

# Results
res_herbivores <- results(DESEQ_herbivores, tidy=TRUE,
                          contrast = c("captive_wild", "captive", "wild"))
View(res_herbivores)

# Table of results for significant ASVs
sigASVs_herbivores <- res_herbivores %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1.5) %>%
  dplyr::rename(ASV = row)
sigASVs_herbivores

# Vector of significant ASVs
sigASVs_herbivores_vec <- sigASVs_herbivores %>%
  pull(ASV)

# Prune phyloseq object
herbivores_filt <- prune_taxa(sigASVs_herbivores_vec, herbivores_samples)

# Add taxonomy onto DESeq results table
merged_results_herbivores <- tax_table(herbivores_filt) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_herbivores, by = "ASV") %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = str_remove(Phylum, "^p__")) %>% # Clean Phylum column
  mutate(Phylum = factor(Phylum, levels = unique(Phylum))) 

# Create DESeq plot
herbivores_DESeq_plot <- ggplot(merged_results_herbivores) +
  geom_bar(aes(x = Phylum, y = log2FoldChange, fill = log2FoldChange > 0), stat = "identity") +
  geom_errorbar(aes(x = Phylum, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.2) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Phylum", y = "Log2 Fold Change (captive/wild)") +
  ggtitle("Differentially Abundant ASVs of Herbivores (Diet Type)") +
  theme(legend.position = "none")
herbivores_DESeq_plot

# Export the phylum plot
ggsave(filename = "DESeq_phylum_plot_herbivores.png", herbivores_DESeq_plot,
       height = 5, width = 15)





############# B) OMNIVORES
# Subset samples to omnivores
omnivore_samples <- subset_samples(zoo_plus1, OchHMC == 'O')

# Subset samples to omnivores
omnivore_plus1 <- transform_sample_counts(omnivore_samples, function(x) x + 1)

# Create DESeq2 object
omnivore_deseq <- phyloseq_to_deseq2(omnivore_plus1, ~ captive_wild)
DESEQ_omnivores <- DESeq(omnivore_deseq)

# Results
res_omnivores <- results(DESEQ_omnivores, tidy=TRUE,
  contrast = c("captive_wild", "captive", "wild"))
View(res_omnivores)

# Table of results for significant ASVs
sigASVs_omnivores <- res_omnivores %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)
sigASVs_omnivores

# Vector of significant ASVs
sigASVs_omnivores_vec <- sigASVs_omnivores %>%
  pull(ASV)

# Prune phyloseq object
omnivore_filt <- prune_taxa(sigASVs_omnivores_vec, omnivore_samples)

# Add taxonomy onto DESeq results table
merged_results_omnivores <- tax_table(omnivore_filt) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_omnivores, by = "ASV") %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = str_remove(Phylum, "^p__")) %>% # Clean Phylum column
  mutate(Phylum = factor(Phylum, levels = unique(Phylum))) 

# Create DESeq plot
omnivores_DESeq_plot <- ggplot(merged_results_omnivores) +
  geom_bar(aes(x = Phylum, y = log2FoldChange, fill = log2FoldChange > 0), stat = "identity") +
  geom_errorbar(aes(x = Phylum, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.2) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Phylum", y = "Log2 Fold Change (captive/wild)") +
  ggtitle("Differentially Abundant ASVs of Omnivores (Diet Type)") +
  theme(legend.position = "none")
omnivores_DESeq_plot

# Export the phylum plot
ggsave(filename = "DESeq_phylum_plot_omnivores.png", omnivores_DESeq_plot,
       height = 5, width = 15)




############# C) CARNIVORES
# Subset samples to carnivores
carnivore_samples <- subset_samples(zoo_plus1, OchHMC == 'C')

# Subset samples to carnivores
carnivore_plus1 <- transform_sample_counts(carnivore_samples, function(x) x + 1)

# Create DESeq2 object
carnivore_deseq <- phyloseq_to_deseq2(carnivore_plus1, ~ captive_wild)
DESEQ_carnivores <- DESeq(carnivore_deseq)

# Results
res_carnivores <- results(DESEQ_carnivores, tidy=TRUE,
                          contrast = c("captive_wild", "captive", "wild"))
View(res_carnivores)

# Table of results for significant ASVs
sigASVs_carnivores <- res_carnivores %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1.5) %>%
  dplyr::rename(ASV = row)
sigASVs_carnivores

# Vector of significant ASVs
sigASVs_carnivores_vec <- sigASVs_carnivores %>%
  pull(ASV)

# Prune phyloseq object
carnivore_filt <- prune_taxa(sigASVs_carnivores_vec, carnivore_samples)

# Add taxonomy onto DESeq results table
merged_results_carnivores <- tax_table(carnivore_filt) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_carnivores, by = "ASV") %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = str_remove(Phylum, "^p__")) %>% # Clean Phylum column
  mutate(Phylum = factor(Phylum, levels = unique(Phylum))) 

# Create DESeq plot
carnivores_DESeq_plot <- ggplot(merged_results_carnivores) +
  geom_bar(aes(x = Phylum, y = log2FoldChange, fill = log2FoldChange > 0), stat = "identity") +
  geom_errorbar(aes(x = Phylum, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.2) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Phylum", y = "Log2 Fold Change (captive/wild)") +
  ggtitle("Differentially Abundant ASVs of Carnivores (Diet Type)") +
  theme(legend.position = "none")
carnivores_DESeq_plot


# Export the phylum plot
ggsave(filename = "DESeq_phylum_plot_carnivores.png", carnivores_DESeq_plot,
       height = 5, width = 15)






#--------------------------- II. Determining the # of points above/below 0

# Add a new column to each merged results table to indicate diet type
merged_results_omnivores <- merged_results_omnivores %>%
  mutate(Diet = "Omnivore")

merged_results_herbivores <- merged_results_herbivores %>%
  mutate(Diet = "Herbivore")

merged_results_carnivores <- merged_results_carnivores %>%
  mutate(Diet = "Carnivore")

# Combine all three data frames into one
combined_results_diet <- bind_rows(merged_results_omnivores, merged_results_herbivores, merged_results_carnivores)

# Calculate mean and standard deviation for each Phylum and Diet type
count_values_diet <- combined_results %>%
  group_by(Diet) %>%
  summarise(
    count_above_0 = sum(log2FoldChange > 0, na.rm = TRUE),
    count_below_0 = sum(log2FoldChange < 0, na.rm = TRUE))
count_values_diet

# Save the merged results as a CSV file
write.csv(count_values_diet, "DESeq_graphvalues_diet.csv", row.names = FALSE)




#--------------------------- III. Determining which phylum most abundant

# Clean Phylum by removing the ".number" part
phylum_abundance_diet <- combined_results %>%
  mutate(Phylum = str_remove(Phylum, "\\..*")) %>%  # Remove ".number" part
  group_by(Diet, Phylum) %>%
  summarise(
    count_above_0 = sum(log2FoldChange > 0, na.rm = TRUE),
    count_below_0 = sum(log2FoldChange < 0, na.rm = TRUE))
phylum_abundance_diet

# Save the merged results as a CSV file
write.csv(phylum_abundance_diet, "DESeq_abundvalues_diet.csv", row.names = FALSE)

