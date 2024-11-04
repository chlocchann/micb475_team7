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

############# A) FOREGUT
# Subset samples to foregut fermenters
foregut_samples <- subset_samples(zoo_plus1, GutFermentersHMC == 'Foregut')

#### DESeq for Captive vs Wild for Foregut Fermenters ####
# Create DESeq2 object
deseq_foregut <- phyloseq_to_deseq2(foregut_samples, ~ captive_wild)
DESEQ_foregut <- DESeq(deseq_foregut)

# Results
res_foregut <- results(DESEQ_foregut, tidy=TRUE, 
                       contrast = c("captive_wild", "captive", "wild"))
View(res_foregut)

# Table of results for significant ASVs
sigASVs_foregut <- res_foregut |>
  filter(padj < 0.05 & abs(log2FoldChange) > 1.5) |>
  dplyr::rename(ASV = row)
View(sigASVs_foregut)

# Vector of significant ASVs
sigASVs_foregut_vec <- sigASVs_foregut %>%
  pull(ASV)

# Identifying the differentially abundant ASVs at the phylum level
foregut_DESeq <- prune_taxa(sigASVs_foregut_vec, foregut_samples)
foregut_sigASVs <- tax_table(foregut_DESeq) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_foregut, by = "ASV") %>%
  arrange(log2FoldChange) %>%
  # Remove 'p__' prefix and any numbers after the phylum name
  mutate(Phylum_Clean = gsub("p__|\\..*", "", Phylum)) %>%
  mutate(Phylum_Clean = factor(Phylum_Clean, levels = unique(Phylum_Clean)))

# Plot at the phylum level
foregut_DESeq_plot <- ggplot(foregut_sigASVs) +
  geom_point(aes(x = Phylum_Clean, y = log2FoldChange)) +
  geom_errorbar(aes(x = Phylum_Clean, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Phylum", y = "Log2 Fold Change (captive/wild)") +
  ggtitle("Foregut (Gut Fermentation)")
foregut_DESeq_plot

# Export the phylum plot
ggsave(filename = "DESeq_phylum_plot_foregut_10000.png", foregut_DESeq_plot,
       height = 5, width = 7)



############# -----------------------------------------------------------


############# B) HINDGUT
hindgut_samples <- subset_samples(zoo_plus1, GutFermentersHMC == 'Hindgut')

#### DESeq for Captive vs Wild for Hindgut Fermenters ####
# Create DESeq2 object
deseq_hindgut <- phyloseq_to_deseq2(hindgut_samples, ~ captive_wild)
DESEQ_hindgut <- DESeq(deseq_hindgut)

# Results
res_hindgut <- results(DESEQ_hindgut, tidy=TRUE, 
                       contrast = c("captive_wild", "captive", "wild"))
View(res_hindgut)

# Table of results
sigASVs_hindgut <- res_hindgut |>
  filter(padj < 0.05 & abs(log2FoldChange) > 1.5) |>
  dplyr::rename(ASV = row)
View(sigASVs_hindgut) #ensure we do have groups with significant difference btw captive & wild

# Vector of significant ASVs
sigASVs_hindgut_vec <- sigASVs_hindgut %>%
  pull(ASV)

# Identifying the differentially abundant ASVs at the phylum level
hindgut_DESeq <- prune_taxa(sigASVs_hindgut_vec, hindgut_samples)

hindgut_sigASVs <- tax_table(hindgut_DESeq) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_hindgut, by = "ASV") %>%
  arrange(log2FoldChange) %>%
  # Remove 'p__' prefix and any numbers after the phylum name
  mutate(Phylum_Clean = gsub("p__|\\..*", "", Phylum)) %>%
  mutate(Phylum_Clean = factor(Phylum_Clean, levels = unique(Phylum_Clean)))

# Plot at the phylum level
hindgut_DESeq_plot <- ggplot(hindgut_sigASVs) +
  geom_point(aes(x = Phylum_Clean, y = log2FoldChange)) +
  geom_errorbar(aes(x = Phylum_Clean, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Phylum", y = "Log2 Fold Change (captive/wild)") +
  ggtitle("Hindgut (Gut Fermentation)")
hindgut_DESeq_plot

# Export the phylum plot
ggsave(filename = "DESeq_phylum_plot_hindgut_10000.png",hindgut_DESeq_plot,
       height = 5, width = 7)




############# -----------------------------------------------------------

############# C) NON-FERM
# Subset samples to non-fermenters
nonfermentor_samples <- subset_samples(zoo_plus1, GutFermentersHMC == 'N')

#### DESeq for Captive vs Wild for Non-Fermentors ####
# Create DESeq2 object
deseq_nonfermentor <- phyloseq_to_deseq2(nonfermentor_samples, ~ captive_wild)
DESEQ_nonfermentor <- DESeq(deseq_nonfermentor)

# Results
res_nonfermentor <- results(DESEQ_nonfermentor, tidy=TRUE, 
                            contrast = c("captive_wild", "captive", "wild"))
View(res_nonfermentor)

# Table of results for significant ASVs
sigASVs_nonfermentor <- res_nonfermentor |>
  filter(padj < 0.05 & abs(log2FoldChange) > 1.5) |>
  dplyr::rename(ASV = row)
View(sigASVs_nonfermentor) #ensure we do have groups with significant difference btw captive & wild

# Vector of significant ASVs
sigASVs_nonfermentor_vec <- sigASVs_nonfermentor %>%
  pull(ASV)

# Identifying the differentially abundant ASVs at the phylum level
nonfermentor_DESeq <- prune_taxa(sigASVs_nonfermentor_vec, nonfermentor_samples)
nonfermentor_sigASVs <- tax_table(nonfermentor_DESeq) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV") %>%
  right_join(sigASVs_nonfermentor, by = "ASV") %>%
  arrange(log2FoldChange) %>%
  # Remove 'p__' prefix and any numbers after the phylum name
  mutate(Phylum_Clean = gsub("p__|\\..*", "", Phylum)) %>%
  mutate(Phylum_Clean = factor(Phylum_Clean, levels = unique(Phylum_Clean)))

# Plot at the phylum level
nonfermentor_DESeq_plot <- ggplot(nonfermentor_sigASVs) +
  geom_point(aes(x = Phylum_Clean, y = log2FoldChange)) +
  geom_errorbar(aes(x = Phylum_Clean, ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE)) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(x = "Phylum", y = "Log2 Fold Change (captive/wild)") +
  ggtitle("Non-Fermentors (Gut Fermentation)")
nonfermentor_DESeq_plot

# Export the phylum plot
ggsave(filename = "DESeq_phylum_plot_nonfermentor_10000.png", nonfermentor_DESeq_plot,
       height = 5, width = 7)







################### II. "POSITIVE EVENTS": noticed some above 0 ASVs for log2FoldChange for categories -------------
# meaning higher in adaptive than wild
# which are relatively (unexpected results)



############# A) FOREGUT
posit_foregut_ASVs <- foregut_sigASVs %>%
  filter(log2FoldChange > 0) %>%
  pull(ASV)

posit_foregut_ASVs
    # [1] "f7f334005597d6d2887e5a1671930134" "1be99e545d159b8da84129d7b2db1d77" "531935cdbf37053fef99397d5cfaa25a"
    # [4] "e9a98931ac582b68717a1d71c4964f3f" "f38b1cb0b64fd29cd983733ecdc2ebe1" "1b938401ec45f73e46f52e7322a7bd47"
    # [7] "cef00af42c4eed7cc366a68d7faa4c27" "79db7e334bfce51b12716bac4d9628a4" "1b00280ce9cd074cd6204adb083fb021"
    # [10] "3dda67b351af428a509efa98f4acda4e" "1748e5b7b634fb98675986bbf10b66d6" "43f3e88dcb2fbce625cf011424a2d01a"
    # [13] "37afc484bc43bc072162eeb576c474f9" "bdbefdad98fb1a9d5b63bea068b4db17"

# Prune the phyloseq object to keep only the specified ASVs
target_ASV_foregut <- c("f7f334005597d6d2887e5a1671930134", "1be99e545d159b8da84129d7b2db1d77", 
                        "531935cdbf37053fef99397d5cfaa25a", "e9a98931ac582b68717a1d71c4964f3f", 
                        "f38b1cb0b64fd29cd983733ecdc2ebe1", "1b938401ec45f73e46f52e7322a7bd47", 
                        "cef00af42c4eed7cc366a68d7faa4c27", "79db7e334bfce51b12716bac4d9628a4",
                        "1b00280ce9cd074cd6204adb083fb021", "3dda67b351af428a509efa98f4acda4e", 
                        "1748e5b7b634fb98675986bbf10b66d6", "43f3e88dcb2fbce625cf011424a2d01a", 
                        "37afc484bc43bc072162eeb576c474f9", "bdbefdad98fb1a9d5b63bea068b4db17")

target_ASV_phyloseq_foregut <- prune_taxa(target_ASV_foregut, foregut_samples)

# Get the taxonomic information for the specified ASVs
tax_info_foregut <- tax_table(target_ASV_phyloseq_foregut) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")
tax_info_foregut

# Get the OTU table to see the counts for the specified ASVs
otu_counts_foregut <- otu_table(target_ASV_phyloseq_foregut) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")
otu_counts_foregut

# Combine taxonomic information with OTU counts
combined_info_foregut <- merge(tax_info_foregut, otu_counts_foregut, by = "ASV")
combined_info_foregut

# Identify the columns for sample counts (assuming they start after the Species column)
count_columns_foregut <- names(combined_info_foregut)[-(1:8)]  # Adjust based on the number of taxonomy columns

# Reshape combined_info to have sample IDs as a column
combined_info_long_foregut <- combined_info_foregut %>%
  pivot_longer(cols = all_of(count_columns_foregut), 
               names_to = "SampleID",  
               values_to = "Count") 

# Get sample data to understand the context of these ASVs
sample_info_foregut <- sample_data(target_ASV_phyloseq_foregut) %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID")

# Merge the reshaped combined_info with sample information and view
combined_results_foregut <- merge(combined_info_long_foregut, sample_info_foregut, by = "SampleID", all.x = TRUE)
View(combined_results_foregut)

# Save the merged results as a CSV file
write.csv(combined_results_foregut, "positive_ASV_foregut.csv", row.names = FALSE)





############# -----------------------------------------------------------

############# B) HINDGUT
posit_hindgut_ASVs <- hindgut_sigASVs %>%
  filter(log2FoldChange > 0) %>%
  pull(ASV)

posit_hindgut_ASVs
    # [1] "60f1d2183c3663fc8d385631d75616ff" "bc511f2c8807bf72a5dd13d5a8dfac1d"

# Prune the phyloseq object to keep only the specified ASVs
target_ASV_hindgut <- c("60f1d2183c3663fc8d385631d75616ff", "bc511f2c8807bf72a5dd13d5a8dfac1d")

target_ASV_phyloseq_hindgut <- prune_taxa(target_ASV_hindgut, hindgut_samples)

# Get the taxonomic information for the specified ASVs
tax_info_hindgut <- tax_table(target_ASV_phyloseq_hindgut) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")
tax_info_hindgut

# Get the OTU table to see the counts for the specified ASVs
otu_counts_hindgut <- otu_table(target_ASV_phyloseq_hindgut) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")
otu_counts_hindgut

# Combine taxonomic information with OTU counts
combined_info_hindgut <- merge(tax_info_hindgut, otu_counts_hindgut, by = "ASV")
combined_info_hindgut

# Identify the columns for sample counts (assuming they start after the Species column)
count_columns_hindgut <- names(combined_info_hindgut)[-(1:8)]  # Adjust based on the number of taxonomy columns

# Reshape combined_info to have sample IDs as a column
combined_info_long_hindgut <- combined_info_hindgut %>%
  pivot_longer(cols = all_of(count_columns_hindgut), 
               names_to = "SampleID",  
               values_to = "Count")

# Get sample data to understand the context of these ASVs
sample_info_hindgut <- sample_data(target_ASV_phyloseq_hindgut) %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID")

# Merge the reshaped combined_info with sample information and view
combined_results_hindgut <- merge(combined_info_long_hindgut, sample_info_hindgut, by = "SampleID", all.x = TRUE)
View(combined_results_hindgut)

# Save the merged results as a CSV file
write.csv(combined_results_hindgut, "positive_ASV_hindgut.csv", row.names = FALSE)



############# -----------------------------------------------------------

############# C) NONFERMENTOR
posit_nonfermentor_ASVs <- nonfermentor_sigASVs %>%
  filter(log2FoldChange > 0) %>%
  pull(ASV)

posit_nonfermentor_ASVs
    # [1] "55b5ad922ce2b83892e0449eaee0dd8f" "055aec5edbd03daea137c9fe44ec23af" "cc86fc399e282b8bfda8ca18779fffac"
    # [4] "e62877c9434952a1d6e72fd6e7cb3b23" "26eb30cfae17a54f95d5b519d330f3aa" "844f46a02a277584eae636be55ac3a91"
    # [7] "6c290419ad5396c5c3f89d18335f06ff" "3b170891f72b36e22a5b615800159838" "ba15b68e1af31016806b38f0898ed018"
    # [10] "c3a7c866f18394e5fc285f4e78bf734b" "3443833956efe96c2c69cc0d0f761f8f" "bf85fc4661d163e33ec8d844c01515ec"
    # [13] "b13d7b642a69f3bc3b8d7f4e6aba0205" "e615ea722dda4acfc3da4648ed487eae" "5db5f38471441cc814e8d20104339fc4"
    # [16] "f001e292e9e19ce85465aaf6413bde55" "b545600dffab681ebe77a5024816a9b5" "a6cd3704b9d621add77f9305332721eb"
    # [19] "d5366f4e78d8e11e5e3ca2134d48f3bd" "be03aa0f984a3378fb41b82da6d69e13" "ca3a12f5eaccca9aa25853a1ab844f63"
    # [22] "7cf74cb58f6a831d55bdc8ffaef348ac" "5afef9c3934a3fd690f59a2675b17ad2" "40420155a410c1d17a505731cd700146"
    # [25] "53d42797b8ac1fa460a02c8a2904581e" "f50defe8cf4a491bc5377919cec7e331"

# Prune the phyloseq object to keep only the specified ASVs
target_ASV_nonfermentor <- c("55b5ad922ce2b83892e0449eaee0dd8f", "055aec5edbd03daea137c9fe44ec23af",
                             "cc86fc399e282b8bfda8ca18779fffac", "e62877c9434952a1d6e72fd6e7cb3b23", 
                             "26eb30cfae17a54f95d5b519d330f3aa", "844f46a02a277584eae636be55ac3a91", 
                             "6c290419ad5396c5c3f89d18335f06ff", "3b170891f72b36e22a5b615800159838", 
                             "ba15b68e1af31016806b38f0898ed018", "c3a7c866f18394e5fc285f4e78bf734b", 
                             "3443833956efe96c2c69cc0d0f761f8f", "bf85fc4661d163e33ec8d844c01515ec", 
                             "b13d7b642a69f3bc3b8d7f4e6aba0205", "e615ea722dda4acfc3da4648ed487eae", 
                             "5db5f38471441cc814e8d20104339fc4", "f001e292e9e19ce85465aaf6413bde55", 
                             "b545600dffab681ebe77a5024816a9b5", "a6cd3704b9d621add77f9305332721eb", 
                             "d5366f4e78d8e11e5e3ca2134d48f3bd", "be03aa0f984a3378fb41b82da6d69e13", 
                             "ca3a12f5eaccca9aa25853a1ab844f63", "7cf74cb58f6a831d55bdc8ffaef348ac", 
                             "5afef9c3934a3fd690f59a2675b17ad2", "40420155a410c1d17a505731cd700146", 
                             "53d42797b8ac1fa460a02c8a2904581e", "f50defe8cf4a491bc5377919cec7e331")

target_ASV_phyloseq_nonfermentor <- prune_taxa(target_ASV_nonfermentor, nonfermentor_samples)

# Get the taxonomic information for the specified ASVs
tax_info_nonfermentor <- tax_table(target_ASV_phyloseq_nonfermentor) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")
tax_info_nonfermentor

# Get the OTU table to see the counts for the specified ASVs
otu_counts_nonfermentor <- otu_table(target_ASV_phyloseq_nonfermentor) %>%
  as.data.frame() %>%
  rownames_to_column(var = "ASV")
otu_counts_nonfermentor

# Combine taxonomic information with OTU counts
combined_info_nonfermentor <- merge(tax_info_nonfermentor, otu_counts_nonfermentor, by = "ASV")
combined_info_nonfermentor

# Identify the columns for sample counts (assuming they start after the Species column)
count_columns_nonfermentor <- names(combined_info_nonfermentor)[-(1:8)]  # Adjust based on the number of taxonomy columns

# Reshape combined_info to have sample IDs as a column
combined_info_long_nonfermentor <- combined_info_nonfermentor %>%
  pivot_longer(cols = all_of(count_columns_nonfermentor), 
               names_to = "SampleID",  
               values_to = "Count")

# Get sample data to understand the context of these ASVs
sample_info_nonfermentor <- sample_data(target_ASV_phyloseq_nonfermentor) %>%
  as.data.frame() %>%
  rownames_to_column(var = "SampleID")

# Merge the reshaped combined_info with sample information and view
combined_results_nonfermentor <- merge(combined_info_long_nonfermentor, sample_info_nonfermentor, by = "SampleID", all.x = TRUE)
View(combined_results_nonfermentor)

# Save the merged results as a CSV file
write.csv(combined_results_nonfermentor, "positive_ASV_nonfermentor.csv", row.names = FALSE)
