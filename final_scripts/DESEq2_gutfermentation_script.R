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

############# A) FOREGUT
# Subset samples to foregut
foregut_samples <- subset_samples(zoo_plus1, GutFermentersHMC == 'Foregut')

# Subset samples to foregut
foregut_plus1 <- transform_sample_counts(foregut_samples, function(x) x + 1)

# Create DESeq2 object
foregut_deseq <- phyloseq_to_deseq2(foregut_plus1, ~ captive_wild)
DESEQ_foregut <- DESeq(foregut_deseq)

# Results
res_foregut <- results(DESEQ_foregut, tidy=TRUE,
                       contrast = c("captive_wild", "captive", "wild"))
View(res_foregut)

# Table of results for significant ASVs
sigASVs_foregut <- res_foregut %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1.5) %>%
  dplyr::rename(ASV = row)
sigASVs_foregut

# Vector of significant ASVs
sigASVs_foregut_vec <- sigASVs_foregut %>%
  pull(ASV)

# Prune phyloseq object
foregut_filt <- prune_taxa(sigASVs_foregut_vec, foregut_samples)

# Add taxonomy onto DESeq results table
merged_results_foregut <- tax_table(foregut_filt) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_foregut, by = "ASV") %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = str_remove(Phylum, "^p__")) %>% # Clean Phylum column
  mutate(Phylum = factor(Phylum, levels = unique(Phylum))) 

# Create DESeq plot
foregut_DESeq_plot <- ggplot(merged_results_foregut) +
  geom_bar(aes(x = Phylum, y = log2FoldChange, fill = log2FoldChange > 0), stat = "identity") +
  geom_errorbar(aes(x = Phylum, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.2) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Phylum", y = "Log2 Fold Change (captive/wild)") +
  ggtitle("Differentially Abundant ASVs of Foregut (Gut Fermentation Type)") +
  theme(legend.position = "none")
foregut_DESeq_plot

# Export the phylum plot
ggsave(filename = "DESeq_phylum_plot_foregut.png", foregut_DESeq_plot,
       height = 5, width = 15)






############# B) HINDGUT
# Subset samples to hindgut
hindgut_samples <- subset_samples(zoo_plus1, GutFermentersHMC == 'Hindgut')

# Subset samples to hindgut
hindgut_plus1 <- transform_sample_counts(hindgut_samples, function(x) x + 1)

# Create DESeq2 object
hindgut_deseq <- phyloseq_to_deseq2(hindgut_plus1, ~ captive_wild)
DESEQ_hindgut <- DESeq(hindgut_deseq)

# Results
res_hindgut <- results(DESEQ_hindgut, tidy=TRUE,
                       contrast = c("captive_wild", "captive", "wild"))
View(res_hindgut)

# Table of results for significant ASVs
sigASVs_hindgut <- res_hindgut %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)
sigASVs_hindgut

# Vector of significant ASVs
sigASVs_hindgut_vec <- sigASVs_hindgut %>%
  pull(ASV)

# Prune phyloseq object
hindgut_filt <- prune_taxa(sigASVs_hindgut_vec, hindgut_samples)

# Add taxonomy onto DESeq results table
merged_results_hindgut <- tax_table(hindgut_filt) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_hindgut, by = "ASV") %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = str_remove(Phylum, "^p__")) %>% # Clean Phylum column
  mutate(Phylum = factor(Phylum, levels = unique(Phylum))) 

# Create DESeq plot
hindgut_DESeq_plot <- ggplot(merged_results_hindgut) +
  geom_bar(aes(x = Phylum, y = log2FoldChange, fill = log2FoldChange > 0), stat = "identity") +
  geom_errorbar(aes(x = Phylum, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.2) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Phylum", y = "Log2 Fold Change (captive/wild)") +
  ggtitle("Differentially Abundant ASVs of Hindgut (Gut Fermentation Type)") +
  theme(legend.position = "none")
hindgut_DESeq_plot

# Export the phylum plot
ggsave(filename = "DESeq_phylum_plot_hindgut.png", hindgut_DESeq_plot,
       height = 5, width = 15)




############# C) NONFERM
# Subset samples to nonferm
nonferm_samples <- subset_samples(zoo_plus1, GutFermentersHMC == 'N')

# Subset samples to nonferm
nonferm_plus1 <- transform_sample_counts(nonferm_samples, function(x) x + 1)

# Create DESeq2 object
nonferm_deseq <- phyloseq_to_deseq2(nonferm_plus1, ~ captive_wild)
DESEQ_nonferm <- DESeq(nonferm_deseq)

# Results
res_nonferm <- results(DESEQ_nonferm, tidy=TRUE,
                       contrast = c("captive_wild", "captive", "wild"))
View(res_nonferm)

# Table of results for significant ASVs
sigASVs_nonferm <- res_nonferm %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1.5) %>%
  dplyr::rename(ASV = row)
sigASVs_nonferm

# Vector of significant ASVs
sigASVs_nonferm_vec <- sigASVs_nonferm %>%
  pull(ASV)

# Prune phyloseq object
nonferm_filt <- prune_taxa(sigASVs_nonferm_vec, nonferm_samples)

# Add taxonomy onto DESeq results table
merged_results_nonferm <- tax_table(nonferm_filt) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_nonferm, by = "ASV") %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = str_remove(Phylum, "^p__")) %>% # Clean Phylum column
  mutate(Phylum = factor(Phylum, levels = unique(Phylum))) 

# Create DESeq plot
nonferm_DESeq_plot <- ggplot(merged_results_nonferm) +
  geom_bar(aes(x = Phylum, y = log2FoldChange, fill = log2FoldChange > 0), stat = "identity") +
  geom_errorbar(aes(x = Phylum, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.2) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Phylum", y = "Log2 Fold Change (captive/wild)") +
  ggtitle("Differentially Abundant ASVs of Nonferm (Gut Fermentation Type)") +
  theme(legend.position = "none")
nonferm_DESeq_plot

# Export the phylum plot
ggsave(filename = "DESeq_phylum_plot_nonferm.png", nonferm_DESeq_plot,
       height = 5, width = 15)




#--------------------------- II. Determining the # of points above/below 0

# Add a new column to each merged results table to indicate fermentation type
merged_results_foregut <- merged_results_foregut %>%
  mutate(Ferm = "Foregut")

merged_results_hindgut <- merged_results_hindgut %>%
  mutate(Ferm = "Hindgut")

merged_results_nonferm <- merged_results_nonferm %>%
  mutate(Ferm = "Nonfermentor")

# Combine all three data frames into one
combined_results_ferm <- bind_rows(merged_results_foregut, merged_results_hindgut, merged_results_nonferm)

# Calculate mean and standard deviation for each Phylum and Fermentation type
count_values_ferm <- combined_results_ferm %>%
  group_by(Ferm) %>%
  summarise(
    count_above_0 = sum(log2FoldChange > 0, na.rm = TRUE),
    count_below_0 = sum(log2FoldChange < 0, na.rm = TRUE))
count_values_ferm

# Save the merged results as a CSV file
write.csv(count_values_ferm, "DESeq_graphvalues_ferm.csv", row.names = FALSE)




#--------------------------- III. Determining which phylum most abundant

# Clean Phylum by removing the ".number" part
phylum_abundance_ferm <- combined_results_ferm %>%
  mutate(Phylum = str_remove(Phylum, "\\..*")) %>%  # Remove ".number" part
  group_by(Ferm, Phylum) %>%
  summarise(
    count_above_0 = sum(log2FoldChange > 0, na.rm = TRUE),
    count_below_0 = sum(log2FoldChange < 0, na.rm = TRUE))
phylum_abundance_ferm

# Save the merged results as a CSV file
write.csv(phylum_abundance_ferm, "DESeq_abundvalues_ferm.csv", row.names = FALSE)
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

############# A) FOREGUT
# Subset samples to foregut
foregut_samples <- subset_samples(zoo_plus1, GutFermentersHMC == 'Foregut')

# Subset samples to foregut
foregut_plus1 <- transform_sample_counts(foregut_samples, function(x) x + 1)

# Create DESeq2 object
foregut_deseq <- phyloseq_to_deseq2(foregut_plus1, ~ captive_wild)
DESEQ_foregut <- DESeq(foregut_deseq)

# Results
res_foregut <- results(DESEQ_foregut, tidy=TRUE,
                       contrast = c("captive_wild", "captive", "wild"))
View(res_foregut)

# Table of results for significant ASVs
sigASVs_foregut <- res_foregut %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1.5) %>%
  dplyr::rename(ASV = row)
sigASVs_foregut

# Vector of significant ASVs
sigASVs_foregut_vec <- sigASVs_foregut %>%
  pull(ASV)

# Prune phyloseq object
foregut_filt <- prune_taxa(sigASVs_foregut_vec, foregut_samples)

# Add taxonomy onto DESeq results table
merged_results_foregut <- tax_table(foregut_filt) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_foregut, by = "ASV") %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = str_remove(Phylum, "^p__")) %>% # Clean Phylum column
  mutate(Phylum = factor(Phylum, levels = unique(Phylum))) 

# Create DESeq plot
foregut_DESeq_plot <- ggplot(merged_results_foregut) +
  geom_bar(aes(x = Phylum, y = log2FoldChange, fill = log2FoldChange > 0), stat = "identity") +
  geom_errorbar(aes(x = Phylum, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.2) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Phylum", y = "Log2 Fold Change (captive/wild)") +
  ggtitle("Differentially Abundant ASVs of Foregut (Gut Fermentation Type)") +
  theme(legend.position = "none")
foregut_DESeq_plot

# Export the phylum plot
ggsave(filename = "DESeq_phylum_plot_foregut.png", foregut_DESeq_plot,
       height = 5, width = 15)






############# B) HINDGUT
# Subset samples to hindgut
hindgut_samples <- subset_samples(zoo_plus1, GutFermentersHMC == 'Hindgut')

# Subset samples to hindgut
hindgut_plus1 <- transform_sample_counts(hindgut_samples, function(x) x + 1)

# Create DESeq2 object
hindgut_deseq <- phyloseq_to_deseq2(hindgut_plus1, ~ captive_wild)
DESEQ_hindgut <- DESeq(hindgut_deseq)

# Results
res_hindgut <- results(DESEQ_hindgut, tidy=TRUE,
                       contrast = c("captive_wild", "captive", "wild"))
View(res_hindgut)

# Table of results for significant ASVs
sigASVs_hindgut <- res_hindgut %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  dplyr::rename(ASV = row)
sigASVs_hindgut

# Vector of significant ASVs
sigASVs_hindgut_vec <- sigASVs_hindgut %>%
  pull(ASV)

# Prune phyloseq object
hindgut_filt <- prune_taxa(sigASVs_hindgut_vec, hindgut_samples)

# Add taxonomy onto DESeq results table
merged_results_hindgut <- tax_table(hindgut_filt) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_hindgut, by = "ASV") %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = str_remove(Phylum, "^p__")) %>% # Clean Phylum column
  mutate(Phylum = factor(Phylum, levels = unique(Phylum))) 

# Create DESeq plot
hindgut_DESeq_plot <- ggplot(merged_results_hindgut) +
  geom_bar(aes(x = Phylum, y = log2FoldChange, fill = log2FoldChange > 0), stat = "identity") +
  geom_errorbar(aes(x = Phylum, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.2) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Phylum", y = "Log2 Fold Change (captive/wild)") +
  ggtitle("Differentially Abundant ASVs of Hindgut (Gut Fermentation Type)") +
  theme(legend.position = "none")
hindgut_DESeq_plot

# Export the phylum plot
ggsave(filename = "DESeq_phylum_plot_hindgut.png", hindgut_DESeq_plot,
       height = 5, width = 15)




############# C) NONFERM
# Subset samples to nonferm
nonferm_samples <- subset_samples(zoo_plus1, GutFermentersHMC == 'N')

# Subset samples to nonferm
nonferm_plus1 <- transform_sample_counts(nonferm_samples, function(x) x + 1)

# Create DESeq2 object
nonferm_deseq <- phyloseq_to_deseq2(nonferm_plus1, ~ captive_wild)
DESEQ_nonferm <- DESeq(nonferm_deseq)

# Results
res_nonferm <- results(DESEQ_nonferm, tidy=TRUE,
                       contrast = c("captive_wild", "captive", "wild"))
View(res_nonferm)

# Table of results for significant ASVs
sigASVs_nonferm <- res_nonferm %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1.5) %>%
  dplyr::rename(ASV = row)
sigASVs_nonferm

# Vector of significant ASVs
sigASVs_nonferm_vec <- sigASVs_nonferm %>%
  pull(ASV)

# Prune phyloseq object
nonferm_filt <- prune_taxa(sigASVs_nonferm_vec, nonferm_samples)

# Add taxonomy onto DESeq results table
merged_results_nonferm <- tax_table(nonferm_filt) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_nonferm, by = "ASV") %>%
  arrange(log2FoldChange) %>%
  mutate(Phylum = make.unique(Phylum)) %>%
  mutate(Phylum = str_remove(Phylum, "^p__")) %>% # Clean Phylum column
  mutate(Phylum = factor(Phylum, levels = unique(Phylum))) 

# Create DESeq plot
nonferm_DESeq_plot <- ggplot(merged_results_nonferm) +
  geom_bar(aes(x = Phylum, y = log2FoldChange, fill = log2FoldChange > 0), stat = "identity") +
  geom_errorbar(aes(x = Phylum, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.2) +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Phylum", y = "Log2 Fold Change (captive/wild)") +
  ggtitle("Differentially Abundant ASVs of Nonferm (Gut Fermentation Type)") +
  theme(legend.position = "none")
nonferm_DESeq_plot

# Export the phylum plot
ggsave(filename = "DESeq_phylum_plot_nonferm.png", nonferm_DESeq_plot,
       height = 5, width = 15)




#--------------------------- II. Determining the # of points above/below 0

# Add a new column to each merged results table to indicate fermentation type
merged_results_foregut <- merged_results_foregut %>%
  mutate(Ferm = "Foregut")

merged_results_hindgut <- merged_results_hindgut %>%
  mutate(Ferm = "Hindgut")

merged_results_nonferm <- merged_results_nonferm %>%
  mutate(Ferm = "Nonfermentor")

# Combine all three data frames into one
combined_results_ferm <- bind_rows(merged_results_foregut, merged_results_hindgut, merged_results_nonferm)

# Calculate mean and standard deviation for each Phylum and Fermentation type
count_values_ferm <- combined_results_ferm %>%
  group_by(Ferm) %>%
  summarise(
    count_above_0 = sum(log2FoldChange > 0, na.rm = TRUE),
    count_below_0 = sum(log2FoldChange < 0, na.rm = TRUE))
count_values_ferm

# Save the merged results as a CSV file
write.csv(count_values_ferm, "DESeq_graphvalues_ferm.csv", row.names = FALSE)




#--------------------------- III. Determining which phylum most abundant

# Clean Phylum by removing the ".number" part
phylum_abundance_ferm <- combined_results_ferm %>%
  mutate(Phylum = str_remove(Phylum, "\\..*")) %>%  # Remove ".number" part
  group_by(Ferm, Phylum) %>%
  summarise(
    count_above_0 = sum(log2FoldChange > 0, na.rm = TRUE),
    count_below_0 = sum(log2FoldChange < 0, na.rm = TRUE))
phylum_abundance_ferm

# Save the merged results as a CSV file
write.csv(phylum_abundance_ferm, "DESeq_abundvalues_ferm.csv", row.names = FALSE)
