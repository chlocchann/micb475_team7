#!/usr/bin/env Rscript

#### Installing and loading packages #### 
install.packages("BiocManager")

# create list for all packages required for ggpicrust2
pkgs <- c("ALDEx2", "SummarizedExperiment", "Biobase", "devtools", 
          "ComplexHeatmap", "BiocGenerics", "metagenomeSeq", 
          "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2")

# install all packages in 'pkg'
# indicate 'n' when asked if want to update all, some or none
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

# install ggpicrust2
install.packages("ggpicrust2")

# load packages
library(readr)
library(ggpicrust2)
library(tibble)
library(tidyverse)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)
library(dplyr)
library(ggrepel)

#### Load Data ####
# importing PICRUSt2 pathway
abundance_file <- "pathway_abundance.tsv"
abundance_data <- read_delim(abundance_file, delim = "\t", col_names = TRUE, trim_ws = TRUE)
abundance_data  =as.data.frame(abundance_data)

# importing metadata file
metadata <- read_delim("zoo_metadata.tsv")

#### Preparing Tables ####
## DIET - CARNIVORE ## 
# replace all non-carnivore variables in OchHMC column with NA
metadata_diet_carnivore <- metadata
metadata_diet_carnivore$OchHMC[metadata_diet_carnivore$OchHMC %in% c("O - CAPTIVE", "O - WILD", "H - CAPTIVE", "H - WILD")] <- NA

# remove NAs in OchHMC
metadata_diet_carnivore = metadata_diet_carnivore[!is.na(metadata_diet_carnivore$OchHMC),]

# filtering abundance table to only include samples that are in the filtered diet and gut fermentation physiology metadata
carnivore = metadata_diet_carnivore$'#SampleID'
carnivore = append(carnivore, "#OTU ID")
abundance_data_carnivore = abundance_data[, colnames(abundance_data) %in% carnivore]

# removing individuals with no data that cause issues downstream with pathways_daa()
abundance_data_carnivore =  abundance_data_carnivore[, colSums(abundance_data_carnivore != 0) > 0]

# ensuring the rownames for filtered abundance tables are empty (required for functions to run)
rownames(abundance_data_carnivore) = NULL

# verify samples in metadata match samples in abundance_data
carnivore_abundance = rownames(t(abundance_data_carnivore[,-1])) 
metadata_diet_carnivore = metadata_diet_carnivore[metadata_diet_carnivore$`#SampleID` %in% carnivore_abundance,] 

## DIET - OMNIVORE ## 
# replace all non-omnivore variables in OchHMC column with NA
metadata_diet_omnivore <- metadata
metadata_diet_omnivore$OchHMC[metadata_diet_omnivore$OchHMC %in% c("C - CAPTIVE", "C - WILD", "H - CAPTIVE", "H - WILD")] <- NA

# remove NAs in OchHMC
metadata_diet_omnivore = metadata_diet_omnivore[!is.na(metadata_diet_omnivore$OchHMC),]

# filtering abundance table to only include samples that are in the filtered diet and gut fermentation physiology metadata
omnivore = metadata_diet_omnivore$'#SampleID'
omnivore = append(omnivore, "#OTU ID")
abundance_data_omnivore = abundance_data[, colnames(abundance_data) %in% omnivore]

# removing individuals with no data that cause issues downstream with pathways_daa()
abundance_data_omnivore =  abundance_data_omnivore[, colSums(abundance_data_omnivore != 0) > 0]

# ensuring the rownames for filtered abundance tables are empty (required for functions to run)
rownames(abundance_data_omnivore) = NULL

# verify samples in metadata match samples in abundance_data
omnivore_abundance = rownames(t(abundance_data_omnivore[,-1])) 
metadata_diet_omnivore = metadata_diet_omnivore[metadata_diet_omnivore$`#SampleID` %in% omnivore_abundance,] 

## DIET - HERBIVORE ## 
# replace all non-herbivore variables in OchHMC column with NA
metadata_diet_herbivore <- metadata
metadata_diet_herbivore$OchHMC[metadata_diet_herbivore$OchHMC %in% c("C - CAPTIVE", "C - WILD", "O - CAPTIVE", "O - WILD")] <- NA

# remove NAs in OchHMC
metadata_diet_herbivore = metadata_diet_herbivore[!is.na(metadata_diet_herbivore$OchHMC),]

# filtering abundance table to only include samples that are in the filtered diet and gut fermentation physiology metadata
herbivore = metadata_diet_herbivore$'#SampleID'
herbivore = append(herbivore, "#OTU ID")
abundance_data_herbivore = abundance_data[, colnames(abundance_data) %in% herbivore]

# removing individuals with no data that cause issues downstream with pathways_daa()
abundance_data_herbivore =  abundance_data_herbivore[, colSums(abundance_data_herbivore != 0) > 0]

# ensuring the rownames for filtered abundance tables are empty (required for functions to run)
rownames(abundance_data_herbivore) = NULL

# verify samples in metadata match samples in abundance_data
herbivore_abundance = rownames(t(abundance_data_herbivore[,-1])) 
metadata_diet_herbivore = metadata_diet_herbivore[metadata_diet_herbivore$`#SampleID` %in% herbivore_abundance,] 

## GUT FERMENTATION - NON-FERMENTERS ## 
# replace all non-non-fermenter variables in GutFermentersHMC column with NA
metadata_gutferment_n <- metadata
metadata_gutferment_n$GutFermentersHMC[metadata_gutferment_n$GutFermentersHMC %in% c("Hindgut - CAPTIVE", "Hindgut - WILD", "Foregut - CAPTIVE", "Foregut - WILD")] <- NA

# remove NAs in GutFermentersHMC
metadata_gutferment_n = metadata_gutferment_n[!is.na(metadata_gutferment_n$GutFermentersHMC),]

# filtering abundance table to only include samples that are in the filtered diet and gut fermentation physiology metadata
non_fermenter = metadata_gutferment_n$'#SampleID'
non_fermenter = append(non_fermenter, "#OTU ID")
abundance_data_non_ferment = abundance_data[, colnames(abundance_data) %in% non_fermenter]

# removing individuals with no data that cause issues downstream with pathways_daa()
abundance_data_non_ferment =  abundance_data_non_ferment[, colSums(abundance_data_non_ferment != 0) > 0]

# ensuring the rownames for filtered abundance tables are empty (required for functions to run)
rownames(abundance_data_non_ferment) = NULL

# verify samples in metadata match samples in abundance_data
non_fermenter_abundance = rownames(t(abundance_data_non_ferment[,-1])) 
metadata_gutferment_n = metadata_gutferment_n[metadata_gutferment_n$`#SampleID` %in% non_fermenter_abundance,]

## GUT FERMENTATION - HINDGUT ## 
# replace all non-hindgut variables in GutFermentersHMC column with NA
metadata_gutferment_h <- metadata
metadata_gutferment_h$GutFermentersHMC[metadata_gutferment_h$GutFermentersHMC %in% c("N - CAPTIVE", "N - WILD", "Foregut - CAPTIVE", "Foregut - WILD")] <- NA

# remove NAs in GutFermentersHMC
metadata_gutferment_h = metadata_gutferment_h[!is.na(metadata_gutferment_h$GutFermentersHMC),]

# filtering abundance table to only include samples that are in the filtered diet and gut fermentation physiology metadata
hindgut = metadata_gutferment_h$'#SampleID'
hindgut = append(hindgut, "#OTU ID")
abundance_data_hindgut = abundance_data[, colnames(abundance_data) %in% hindgut]

# removing individuals with no data that cause issues downstream with pathways_daa()
abundance_data_hindgut =  abundance_data_hindgut[, colSums(abundance_data_hindgut != 0) > 0]

# ensuring the rownames for filtered abundance tables are empty (required for functions to run)
rownames(abundance_data_hindgut) = NULL

# verify samples in metadata match samples in abundance_data
hindgut_fermenter_abundance = rownames(t(abundance_data_hindgut[,-1])) 
metadata_gutferment_h = metadata_gutferment_h[metadata_gutferment_h$`#SampleID` %in% hindgut_fermenter_abundance,]

## GUT FERMENTATION - FOREGUT ## 
# replace all non-foregut variables in GutFermentersHMC column with NA
metadata_gutferment_f <- metadata
metadata_gutferment_f$GutFermentersHMC[metadata_gutferment_f$GutFermentersHMC %in% c("N - CAPTIVE", "N - WILD", "Hindgut - CAPTIVE", "Hindgut - WILD")] <- NA

# remove NAs in GutFermentersHMC
metadata_gutferment_f = metadata_gutferment_f[!is.na(metadata_gutferment_f$GutFermentersHMC),]

# filtering abundance table to only include samples that are in the filtered diet and gut fermentation physiology metadata
foregut = metadata_gutferment_f$'#SampleID'
foregut = append(foregut, "#OTU ID")
abundance_data_foregut = abundance_data[, colnames(abundance_data) %in% foregut]

# removing individuals with no data that cause issues downstream with pathways_daa()
abundance_data_foregut =  abundance_data_foregut[, colSums(abundance_data_foregut != 0) > 0]

# ensuring the rownames for filtered abundance tables are empty (required for functions to run)
rownames(abundance_data_foregut) = NULL

# verify samples in metadata match samples in abundance_data
foregut_fermenter_abundance = rownames(t(abundance_data_foregut[,-1])) 
metadata_gutferment_f = metadata_gutferment_f[metadata_gutferment_f$`#SampleID` %in% foregut_fermenter_abundance,]

#### DESEq - Heatmaps & PCA Plots ####
## DIET - CARNIVORE ## 
# perform pathway DAA using DESEQ2 method
abundance_daa_carnivore <- pathway_daa(abundance = abundance_data_carnivore %>% column_to_rownames("#OTU ID"), 
                                  metadata = metadata_diet_carnivore, group = "OchHMC", daa_method = "DESeq2")

# annotate MetaCyc pathway so they are more descriptive
metacyc_daa_carnivore <- pathway_annotation(pathway = "MetaCyc", 
                                       daa_results_df = abundance_daa_carnivore, ko_to_kegg = FALSE)

# filter p-values to only significant ones
carnivore_with_p_0.05 <- abundance_daa_carnivore %>% filter(p_values < 0.05)

# changing #OTU ID column to description for the results 
carnivore_desc = inner_join(carnivore_with_p_0.05,metacyc_daa_carnivore, by = "feature")
carnivore_desc$feature = carnivore_desc$description
carnivore_desc = carnivore_desc[,c(1:7)]
colnames(carnivore_desc) = colnames(carnivore_with_p_0.05)

# changing #OTU ID column to description for the abundance table
carnivore_abundance = abundance_data_carnivore %>% filter(`#OTU ID` %in% carnivore_with_p_0.05$feature)
colnames(carnivore_abundance)[1] = "feature"
carnivore_abundance_desc = inner_join(carnivore_abundance,metacyc_daa_carnivore, by = "feature")
carnivore_abundance_desc$feature = carnivore_abundance_desc$description
carnivore_abundance_desc = carnivore_abundance_desc[,-c(20:ncol(carnivore_abundance_desc))] 

# generate a heatmap
CARNIVORE_HEATMAP <- pathway_heatmap(abundance = carnivore_abundance_desc %>% column_to_rownames("feature"), metadata = metadata_diet_carnivore, group = "OchHMC")

# generate pathway PCA plot
CARNIVORE_PCA <- pathway_pca(abundance = abundance_data_carnivore %>% column_to_rownames("#OTU ID"), metadata = metadata_diet_carnivore, group = "OchHMC")

## DIET - OMNIVORE ## 
# perform pathway DAA using DESEQ2 method
abundance_daa_omnivore <- pathway_daa(abundance = abundance_data_omnivore %>% column_to_rownames("#OTU ID"), 
                                       metadata = metadata_diet_omnivore, group = "OchHMC", daa_method = "DESeq2")

# annotate MetaCyc pathway so they are more descriptive
metacyc_daa_omnivore <- pathway_annotation(pathway = "MetaCyc", 
                                            daa_results_df = abundance_daa_omnivore, ko_to_kegg = FALSE)

# filter p-values to only significant ones
omnivore_with_p_0.05 <- abundance_daa_omnivore %>% filter(p_values < 0.05)

# changing #OTU ID column to description for the results 
omnivore_desc = inner_join(omnivore_with_p_0.05,metacyc_daa_omnivore, by = "feature")
omnivore_desc$feature = omnivore_desc$description
omnivore_desc = omnivore_desc[,c(1:7)]
colnames(omnivore_desc) = colnames(omnivore_with_p_0.05)

# changing #OTU ID column to description for the abundance table
omnivore_abundance = abundance_data_omnivore %>% filter(`#OTU ID` %in% omnivore_with_p_0.05$feature)
colnames(omnivore_abundance)[1] = "feature"
omnivore_abundance_desc = inner_join(omnivore_abundance,metacyc_daa_omnivore, by = "feature")
omnivore_abundance_desc$feature = omnivore_abundance_desc$description
omnivore_abundance_desc = omnivore_abundance_desc[,-c(8:ncol(omnivore_abundance_desc))] 

# generate a heatmap
OMNIVORE_HEATMAP <- pathway_heatmap(abundance = omnivore_abundance_desc %>% column_to_rownames("feature"), metadata = metadata_diet_omnivore, group = "OchHMC")
# ERROR: BOTH NO WILD OMNIVORE SAMPLES SHOWN

# generate pathway PCA plot
OMNIVORE_PCA <- pathway_pca(abundance = abundance_data_omnivore %>% column_to_rownames("#OTU ID"), metadata = metadata_diet_omnivore, group = "OchHMC")
# ERROR: ZERO VARIANCE BETWEEN CAPTIVE & WILD OMNIVORES

## DIET - HERBIVORE ## 
# perform pathway DAA using DESEQ2 method
abundance_daa_herbivore <- pathway_daa(abundance = abundance_data_herbivore %>% column_to_rownames("#OTU ID"), 
                                       metadata = metadata_diet_herbivore, group = "OchHMC", daa_method = "DESeq2")

# annotate MetaCyc pathway so they are more descriptive
metacyc_daa_herbivore <- pathway_annotation(pathway = "MetaCyc", 
                                            daa_results_df = abundance_daa_herbivore, ko_to_kegg = FALSE)

# filter p-values to only significant ones
herbivore_with_p_0.05 <- abundance_daa_herbivore %>% filter(p_values < 0.05)

# changing #OTU ID column to description for the results 
herbivore_desc = inner_join(herbivore_with_p_0.05,metacyc_daa_herbivore, by = "feature")
herbivore_desc$feature = herbivore_desc$description
herbivore_desc = herbivore_desc[,c(1:7)]
colnames(herbivore_desc) = colnames(herbivore_with_p_0.05)

# changing #OTU ID column to description for the abundance table
herbivore_abundance = abundance_data_herbivore %>% filter(`#OTU ID` %in% herbivore_with_p_0.05$feature)
colnames(herbivore_abundance)[1] = "feature"
herbivore_abundance_desc = inner_join(herbivore_abundance,metacyc_daa_herbivore, by = "feature")
herbivore_abundance_desc$feature = herbivore_abundance_desc$description
herbivore_abundance_desc = herbivore_abundance_desc[,-c(63:ncol(herbivore_abundance_desc))] 

# generate a heatmap
HERBIVORE_HEATMAP <- pathway_heatmap(abundance = sig_abun_carnivore %>% column_to_rownames("feature"), metadata = metadata_diet_herbivore, group = "OchHMC")

# generate pathway PCA plot
HERBIVORE_PCA <- pathway_pca(abundance = abundance_data_herbivore %>% column_to_rownames("#OTU ID"), metadata = metadata_diet_herbivore, group = "OchHMC")
# ERROR: ZERO VARIANCE BETWEEN CAPTIVE & WILD HERBIVORES

## GUT FERMENTATION - NON-FERMENTERS ## 
# perform pathway DAA using DESEQ2 method
abundance_daa_non_fermenter <- pathway_daa(abundance = abundance_data_non_ferment %>% column_to_rownames("#OTU ID"), 
                                           metadata = metadata_gutferment_n, group = "GutFermentersHMC", daa_method = "DESeq2")

# annotate MetaCyc pathway so they are more descriptive
metacyc_daa_non_fermenter <- pathway_annotation(pathway = "MetaCyc", 
                                                daa_results_df = abundance_daa_non_fermenter, ko_to_kegg = FALSE)

# filter p-values to only significant ones
non_fermenter_with_p_0.05 <- abundance_daa_non_fermenter %>% filter(p_values < 0.05)

# changing #OTU ID column to description for the results 
non_fermenter_desc = inner_join(non_fermenter_with_p_0.05,metacyc_daa_non_fermenter, by = "feature")
non_fermenter_desc$feature = non_fermenter_desc$description
non_fermenter_desc = non_fermenter_desc[,c(1:7)]
colnames(non_fermenter_desc) = colnames(non_fermenter_with_p_0.05)

# changing #OTU ID column to description for the abundance table
non_fermenter_abundance = abundance_data_non_ferment %>% filter(`#OTU ID` %in% non_fermenter_with_p_0.05$feature)
colnames(non_fermenter_abundance)[1] = "feature"
non_fermenter_abundance_desc = inner_join(non_fermenter_abundance,metacyc_daa_non_fermenter, by = "feature")
non_fermenter_abundance_desc$feature = non_fermenter_abundance_desc$description
non_fermenter_abundance_desc = non_fermenter_abundance_desc[,-c(34:ncol(non_fermenter_abundance_desc))] 

# generate a heatmap
NON_FERMENTER_HEATMAP <- pathway_heatmap(abundance = non_fermenter_abundance_desc %>% column_to_rownames("feature"), metadata = metadata_gutferment_n, group = "GutFermentersHMC")

# generate pathway PCA plot
NON_FERMENTER_PCA <- pathway_pca(abundance = abundance_data_non_ferment %>% column_to_rownames("#OTU ID"), metadata = metadata_gutferment_n, group = "GutFermentersHMC")

## GUT FERMENTATION - HINDGUT ## 
# perform pathway DAA using DESEQ2 method
abundance_daa_hindgut <- pathway_daa(abundance = abundance_data_hindgut %>% column_to_rownames("#OTU ID"), 
                                           metadata = metadata_gutferment_h, group = "GutFermentersHMC", daa_method = "DESeq2")

# annotate MetaCyc pathway so they are more descriptive
metacyc_daa_hindgut <- pathway_annotation(pathway = "MetaCyc", 
                                                daa_results_df = abundance_daa_hindgut, ko_to_kegg = FALSE)

# filter p-values to only significant ones
hindgut_with_p_0.05 <- abundance_daa_hindgut %>% filter(p_values < 0.05)

# changing #OTU ID column to description for the results 
hindgut_desc = inner_join(hindgut_with_p_0.05,metacyc_daa_hindgut, by = "feature")
hindgut_desc$feature = hindgut_desc$description
hindgut_desc = hindgut_desc[,c(1:7)]
colnames(hindgut_desc) = colnames(hindgut_with_p_0.05)

# changing #OTU ID column to description for the abundance table
hindgut_abundance = abundance_data_hindgut %>% filter(`#OTU ID` %in% hindgut_with_p_0.05$feature)
colnames(hindgut_abundance)[1] = "feature"
hindgut_abundance_desc = inner_join(hindgut_abundance,metacyc_daa_hindgut, by = "feature")
hindgut_abundance_desc$feature = hindgut_abundance_desc$description
hindgut_abundance_desc = hindgut_abundance_desc[,-c(33:ncol(hindgut_abundance_desc))] 

# generate a heatmap
HINDGUT_HEATMAP <- pathway_heatmap(abundance = hindgut_abundance_desc %>% column_to_rownames("feature"), metadata = metadata_gutferment_h, group = "GutFermentersHMC")

# generate pathway PCA plot
HINDGUT_PCA <- pathway_pca(abundance = abundance_data_hindgut %>% column_to_rownames("#OTU ID"), metadata = metadata_gutferment_h, group = "GutFermentersHMC")
# ERROR: ZERO VARIANCE BETWEEN CAPTIVE & WILD HINDGUT FERMENTERS

## GUT FERMENTATION - FOREGUT ## 
# perform pathway DAA using DESEQ2 method
abundance_daa_foregut <- pathway_daa(abundance = abundance_data_foregut %>% column_to_rownames("#OTU ID"), 
                                  metadata = metadata_gutferment_f, group = "GutFermentersHMC", daa_method = "DESeq2")

# annotate MetaCyc pathway so they are more descriptive
metacyc_daa_foregut <- pathway_annotation(pathway = "MetaCyc", 
                                       daa_results_df = abundance_daa_foregut, ko_to_kegg = FALSE)

# filter p-values to only significant ones
foregut_with_p_0.05 <- abundance_daa_foregut %>% filter(p_values < 0.05)

# changing #OTU ID column to description for the results 
foregut_desc = inner_join(foregut_with_p_0.05,metacyc_daa_foregut, by = "feature")
foregut_desc$feature = foregut_desc$description
foregut_desc = foregut_desc[,c(1:7)]
colnames(foregut_desc) = colnames(foregut_with_p_0.05)

# changing #OTU ID column to description for the abundance table
foregut_abundance = abundance_data_foregut %>% filter(`#OTU ID` %in% foregut_with_p_0.05$feature)
colnames(foregut_abundance)[1] = "feature"
foregut_abundance_desc = inner_join(foregut_abundance,metacyc_daa_foregut, by = "feature")
foregut_abundance_desc$feature = foregut_abundance_desc$description
foregut_abundance_desc = foregut_abundance_desc[,-c(24:ncol(foregut_abundance_desc))] 

# generate a heatmap
FOREGUT_HEATMAP <- pathway_heatmap(abundance = foregut_abundance_desc %>% column_to_rownames("feature"), metadata = metadata_gutferment_f, group = "GutFermentersHMC")

# generate pathway PCA plot
FOREGUT_PCA <- pathway_pca(abundance = abundance_data_foregut %>% column_to_rownames("#OTU ID"), metadata = metadata_gutferment_f, group = "GutFermentersHMC")
# ERROR: ZERO VARIANCE BETWEEN CAPTIVE & WILD FOREGUT FERMENTERS

#### DESEq - Log2FC & Volcano Plots ####
# leading in custom DESEq2 function
source("PICRUSt2_DESeq2_function.R")

## DIET - CARNIVORE ## 
# running custom DESEq2 function
carnivore_res =  DEseq2_function(abundance_data_carnivore, metadata_diet_carnivore, "OchHMC")
carnivore_res$feature =rownames(carnivore_res)
carnivore_res_desc = inner_join(carnivore_res,metacyc_daa_carnivore, by = "feature")
carnivore_res_desc = carnivore_res_desc[, -c(8:13)]

# filter to only include significant pathways based on p-value & log2foldchange
carnivore_sig_res = carnivore_res_desc %>%
  filter(pvalue < 0.05) %>%
  filter(abs(log2FoldChange) > 2) 
carnivore_sig_res <- carnivore_sig_res[order(carnivore_sig_res$log2FoldChange),]

# NOTE: only downregulated pathways found, so filtered for top 10 downregulated pathways
# sorting based on absolute log2FC for top 10 downregulated pathways 
bottom_10_carnivore <- carnivore_sig_res %>%
  filter(log2FoldChange < 0) %>%
  arrange(abs(log2FoldChange)) %>%
  head(10)

# generating a plot of top 10 downregulated pathways
CARNIVORE_LOG2FC <- ggplot(data = bottom_10_carnivore, aes(y = reorder(description, log2FoldChange), x = log2FoldChange, fill = pvalue)) +
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

# creating an additional column to indicate which pathways are upregulated/downregulated/insignificant
carnivore_res_desc <- carnivore_res_desc %>%
  mutate(regulation = case_when(
    pvalue < 0.05 & log2FoldChange > 2  ~ "Upregulated",
    pvalue < 0.05 & log2FoldChange < -2 ~ "Downregulated",
    TRUE ~ "Insignificant"))

# generating a volcano plot showing significant upregulated/downregulated pathways
CARNIVORE_VOLCANO <- ggplot(carnivore_res_desc, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(color = regulation), alpha = 0.7) +  
  labs(x = "Log2 Fold Change",
       y = "-Log10(p-value)",
       color = "Regulatory Status") +  
  theme_minimal() +
  theme(legend.position = "right") +  
  geom_vline(xintercept = c(-2, 2), linetype = "solid", color = "black") + 
  geom_hline(yintercept = -log10(0.05), linetype = "solid", color = "black") +   
  scale_color_manual(values = c("Upregulated" = "indianred",  
                                "Downregulated" = "skyblue",  
                                "Insignificant" = "grey"))  

## DIET - OMNIVORE ## 
# running custom DESEq2 function
omnivore_res =  DEseq2_function(abundance_data_omnivore, metadata_diet_omnivore, "OchHMC")
omnivore_res$feature =rownames(omnivore_res)
omnivore_res_desc = inner_join(omnivore_res,metacyc_daa_omnivore, by = "feature")
omnivore_res_desc = omnivore_res_desc[, -c(8:13)]

# filter to only include significant pathways based on p-value & log2foldchange
omnivore_sig_res = omnivore_res_desc %>%
  filter(pvalue < 0.05) %>%
  filter(abs(log2FoldChange) > 2) 
omnivore_sig_res <- omnivore_sig_res[order(omnivore_sig_res$log2FoldChange),]

# NOTE: only downregulated pathways found, so filtered for top 10 downregulated pathways
# sorting based on absolute log2FC for bottom 10 pathways 
bottom_10_omnivore <- omnivore_sig_res %>%
  filter(log2FoldChange < 0) %>%
  arrange(abs(log2FoldChange)) %>%
  head(10)

# generating a plot of top 10 downregulated pathways
OMNIVORE_LOG2FC <- ggplot(data = bottom_10_omnivore, aes(y = reorder(description, log2FoldChange), x = log2FoldChange, fill = pvalue)) +
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

# creating an additional column to indicate which pathways are upregulated/downregulated/insignificant
omnivore_res_desc <- omnivore_res_desc %>%
  mutate(regulation = case_when(
    pvalue < 0.05 & log2FoldChange > 2  ~ "Upregulated",
    pvalue < 0.05 & log2FoldChange < -2 ~ "Downregulated",
    TRUE ~ "Insignificant"))

# generating a volcano plot showing significant upregulated/downregulated pathways
OMNIVORE_VOLCANO <- ggplot(omnivore_res_desc, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(color = regulation), alpha = 0.7) +  
  labs(x = "Log2 Fold Change",
       y = "-Log10(p-value)",
       color = "Regulatory Status") +  
  theme_minimal() +
  theme(legend.position = "right") +  
  geom_vline(xintercept = c(-2, 2), linetype = "solid", color = "black") + 
  geom_hline(yintercept = -log10(0.05), linetype = "solid", color = "black") +   
  scale_color_manual(values = c("Upregulated" = "indianred",  
                                "Downregulated" = "skyblue",  
                                "Insignificant" = "grey"))  

## DIET - HERBIVORE ## 
# running custom DESEq2 function
herbivore_res =  DEseq2_function(abundance_data_herbivore, metadata_diet_herbivore, "OchHMC")
herbivore_res$feature =rownames(herbivore_res)
herbivore_res_desc = inner_join(herbivore_res,metacyc_daa_herbivore, by = "feature")
herbivore_res_desc = herbivore_res_desc[, -c(8:13)]

# filter to only include significant pathways based on p-value & log2foldchange
herbivore_sig_res = herbivore_res_desc %>%
  filter(pvalue < 0.05) %>%
  filter(abs(log2FoldChange) > 2) 
herbivore_sig_res <- herbivore_sig_res[order(herbivore_sig_res$log2FoldChange),]

# NOTE: only 3 pathways were up/downregulated so no further filtering was performed
# generating a plot of top upregulated & downregulated genes
HERBIVORE_LOG2FC <- ggplot(data = herbivore_sig_res, aes(y = reorder(description, log2FoldChange), x = log2FoldChange, fill = pvalue)) +
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

# creating an additional column to indicate which pathways are upregulated/downregulated/insignificant
herbivore_res_desc <- herbivore_res_desc %>%
  mutate(regulation = case_when(
    pvalue < 0.05 & log2FoldChange > 2  ~ "Upregulated",
    pvalue < 0.05 & log2FoldChange < -2 ~ "Downregulated",
    TRUE ~ "Insignificant"))

# generating a volcano plot showing significant upregulated/downregulated pathways
HERBIVORE_VOLCANO <- ggplot(herbivore_res_desc, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(color = regulation), alpha = 0.7) +  
  labs(x = "Log2 Fold Change",
       y = "-Log10(p-value)",
       color = "Regulatory Status") +  
  theme_minimal() +
  theme(legend.position = "right") +  
  geom_vline(xintercept = c(-2, 2), linetype = "solid", color = "black") + 
  geom_hline(yintercept = -log10(0.05), linetype = "solid", color = "black") +   
  scale_color_manual(values = c("Upregulated" = "indianred",  
                                "Downregulated" = "skyblue",  
                                "Insignificant" = "grey"))  

## GUT FERMENTATION - NON-FERMENTER ## 
# running custom DESEq2 function
nonferment_res =  DEseq2_function(abundance_data_non_ferment, metadata_gutferment_n, "GutFermentersHMC")
nonferment_res$feature =rownames(nonferment_res)
nonferment_res_desc = inner_join(nonferment_res,metacyc_daa_non_fermenter, by = "feature")
nonferment_res_desc = nonferment_res_desc[, -c(8:13)]

# filter to only include significant pathways based on p-value & log2foldchange
nonferment_sig_res = nonferment_res_desc %>%
  filter(pvalue < 0.05) %>%
  filter(abs(log2FoldChange) > 2) 
nonferment_sig_res <- nonferment_sig_res[order(nonferment_sig_res$log2FoldChange),]

# sorting based on absolute log2FC for top upregulated pathways
top_10_nonferment <- nonferment_sig_res %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(10)

top_10_nonferment <- top_10_nonferment %>%
  mutate(group = "Top Upregulated")

# sorting based on absolute log2FC for top downregulated pathways 
bottom_10_nonferment <- nonferment_sig_res %>%
  filter(log2FoldChange < 0) %>%
  arrange(abs(log2FoldChange)) %>%
  head(10)

bottom_10_nonferment <- bottom_10_nonferment %>%
  mutate(group = "Bottom Downregulated")

# combining the top and bottom 10 pathways to plot
sig20_nonferment <- bind_rows(top_10_nonferment, bottom_10_nonferment) %>%
  arrange(abs(log2FoldChange))

# generating a plot of top 10 upregulated & downregulated genes
NON_FERMENTER_LOG2FC <- ggplot(data = sig20_nonferment, aes(y = reorder(description, log2FoldChange), x = log2FoldChange, fill = pvalue)) +
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

# creating an additional column to indicate which pathways are upregulated/downregulated/insignificant
nonferment_res_desc <- nonferment_res_desc %>%
  mutate(regulation = case_when(
    pvalue < 0.05 & log2FoldChange > 2  ~ "Upregulated",
    pvalue < 0.05 & log2FoldChange < -2 ~ "Downregulated",
    TRUE ~ "Insignificant"))

# generating a volcano plot showing significant upregulated/downregulated pathways
NON_FERMENTER_VOLCANO <- ggplot(nonferment_res_desc, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(color = regulation), alpha = 0.7) +  
  labs(x = "Log2 Fold Change",
       y = "-Log10(p-value)",
       color = "Regulatory Status") +  
  theme_minimal() +
  theme(legend.position = "right") +  
  geom_vline(xintercept = c(-2, 2), linetype = "solid", color = "black") + 
  geom_hline(yintercept = -log10(0.05), linetype = "solid", color = "black") +   
  scale_color_manual(values = c("Upregulated" = "indianred",  
                                "Downregulated" = "skyblue",  
                                "Insignificant" = "grey")) 

## GUT FERMENTATION - HINDGUT ## 
# running custom DESEq2 function
hindgut_res =  DEseq2_function(abundance_data_hindgut, metadata_gutferment_h, "GutFermentersHMC")
hindgut_res$feature =rownames(hindgut_res)
hindgut_res_desc = inner_join(hindgut_res,metacyc_daa_hindgut, by = "feature")
hindgut_res_desc = hindgut_res_desc[, -c(8:13)]

# filter to only include significant pathways based on p-value & log2foldchange
hindgut_sig_res = hindgut_res_desc %>%
  filter(pvalue < 0.05) %>%
  filter(abs(log2FoldChange) > 2) 
hindgut_sig_res <- hindgut_sig_res[order(hindgut_sig_res$log2FoldChange),]

# sorting based on absolute log2FC for top 10 pathways 
top_10_hindgut <- hindgut_sig_res %>%
  filter(log2FoldChange > 0) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(10)

top_10_hindgut <- top_10_hindgut %>%
  mutate(group = "Top Upregulated")

# sorting based on absolute log2FC for bottom 10 pathways 
bottom_10_hindgut <- hindgut_sig_res %>%
  filter(log2FoldChange < 0) %>%
  arrange(abs(log2FoldChange)) %>%
  head(10)

bottom_10_hindgut <- bottom_10_hindgut %>%
  mutate(group = "Top Downregulated")

# combining the top and bottom 10 pathways to plot
sig20_hindgut <- bind_rows(top_10_hindgut, bottom_10_hindgut)

# generating a plot of top 10 upregulated & downregulated genes
HINDGUT_LOG2FC <- ggplot(data = sig20_hindgut, aes(y = reorder(description, log2FoldChange), x = log2FoldChange, fill = pvalue)) +
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

# creating an additional column to indicate which pathways are upregulated/downregulated/insignificant
hindgut_res_desc <- hindgut_res_desc %>%
  mutate(regulation = case_when(
    pvalue < 0.05 & log2FoldChange > 2  ~ "Upregulated",
    pvalue < 0.05 & log2FoldChange < -2 ~ "Downregulated",
    TRUE ~ "Insignificant"))

# generating a volcano plot showing significant upregulated/downregulated pathways
HINDGUT_VOLCANO <- ggplot(hindgut_res_desc, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(color = regulation), alpha = 0.7) +  
  labs(x = "Log2 Fold Change",
       y = "-Log10(p-value)",
       color = "Regulatory Status") +  
  theme_minimal() +
  theme(legend.position = "right") +  
  geom_vline(xintercept = c(-2, 2), linetype = "solid", color = "black") + 
  geom_hline(yintercept = -log10(0.05), linetype = "solid", color = "black") +   
  scale_color_manual(values = c("Upregulated" = "indianred",  
                                "Downregulated" = "skyblue",  
                                "Insignificant" = "grey")) 

## GUT FERMENTATION - FOREGUT ## 
# running custom DESEq2 function
foregut_res =  DEseq2_function(abundance_data_foregut, metadata_gutferment_f, "GutFermentersHMC")
foregut_res$feature =rownames(foregut_res)
foregut_res_desc = inner_join(foregut_res,metacyc_daa_foregut, by = "feature")
foregut_res_desc = foregut_res_desc[, -c(8:13)]

# filter to only include significant pathways based on p-value & log2foldchange
foregut_sig_res = foregut_res_desc %>%
  filter(pvalue < 0.05) %>%
  filter(abs(log2FoldChange) > 2) 
foregut_sig_res <- foregut_sig_res[order(foregut_sig_res$log2FoldChange),]

# NOTE: only 3 downregulated pathways were identified so no further filtering was performed
# generating a plot of top downregulated pathways
FOREGUT_LOG2FC <- ggplot(data = foregut_sig_res, aes(y = reorder(description, sort(as.numeric(log2FoldChange))), x= log2FoldChange, fill = pvalue))+
  geom_bar(stat = "identity")+ 
  theme_bw()+
  labs(x = "Log2FoldChange", y="Pathways")

# creating an additional column to indicate which pathways are upregulated/downregulated/insignificant
foregut_res_desc <- foregut_res_desc %>%
  mutate(regulation = case_when(
    pvalue < 0.05 & log2FoldChange > 2  ~ "Upregulated",
    pvalue < 0.05 & log2FoldChange < -2 ~ "Downregulated",
    TRUE ~ "Insignificant"))

# generating a volcano plot showing significant upregulated/downregulated pathways
FOREGUT_VOLCANO <- ggplot(foregut_res_desc, aes(x = log2FoldChange, y = -log10(pvalue))) + 
  geom_point(aes(color = regulation), alpha = 0.7) +  
  labs(x = "Log2 Fold Change",
       y = "-Log10(p-value)",
       color = "Regulatory Status") +  
  theme_minimal() +
  theme(legend.position = "right") +  
  geom_vline(xintercept = c(-2, 2), linetype = "solid", color = "black") + 
  geom_hline(yintercept = -log10(0.05), linetype = "solid", color = "black") +  
  scale_color_manual(values = c("Upregulated" = "indianred",  
                                "Downregulated" = "skyblue",  
                                "Insignificant" = "grey")) 

#### Saving Plots ####
## DIET - CARNIVORE ## 
# saving heatmap
ggsave("CARNIVORE_HEATMAP.PNG"
       , CARNIVORE_HEATMAP
       , height=15, width =12)
# saving pca plot
ggsave("CARNIVORE_PCA.PNG"
       , CARNIVORE_PCA
       , height=6, width =12)
# saving top 10 upregulated & downregulated pathways
ggsave("CARNIVORE_LOG2FC.PNG"
       , CARNIVORE_LOG2FC
       , height=6, width =12)
# saving volcano plot
ggsave("CARNIVORE_VOLCANO.PNG"
       , CARNIVORE_VOLCANO
       , height=6, width =12)

## DIET - OMNIVORE ## 
# saving heatmap
ggsave("OMNIVORE_HEATMAP.PNG"
       , OMNIVORE_HEATMAP
       , height=6, width =12)
# no pca plot to save due to lack of variance between captive vs. wild samples
# saving top 10 upregulated & downregulated pathways
ggsave("OMNIVORE_LOG2FC.PNG"
       , OMNIVORE_LOG2FC
       , height=6, width =12)
# saving volcano plot
ggsave("OMNIVORE_VOLCANO.PNG"
       , OMNIVORE_VOLCANO
       , height=6, width =12)

## DIET - HERBIVORE ## 
# saving heatmap
ggsave("HERBIVORE_HEATMAP.PNG"
       , HERBIVORE_HEATMAP
       , height=6, width =24)
# no pca plot to save due to lack of variance between captive vs. wild samples
# saving top 10 upregulated & downregulated pathways
ggsave("HERBIVORE_LOG2FC.PNG"
       , HERBIVORE_LOG2FC
       , height=6, width =12)
# saving volcano plot
ggsave("HERBIVORE_VOLCANO.PNG"
       , HERBIVORE_VOLCANO
       , height=6, width =12)

## GUT FERMENTATION - NON-FERMENTER ## 
ggsave("NON_FERMENTER_HEATMAP.PNG"
       , NON_FERMENTER_HEATMAP
       , height=15, width =12)
# saving pca plot
ggsave("NON_FERMENTER_PCA.PNG"
       , NON_FERMENTER_PCA
       , height=6, width =12)
# saving top 10 upregulated & downregulated pathways
ggsave("NON_FERMENTER_LOG2FC.PNG"
       , NON_FERMENTER_LOG2FC
       , height=6, width =12)
# saving volcano plot
ggsave("NON_FERMENTER_VOLCANO.PNG"
       , NON_FERMENTER_VOLCANO
       , height=6, width =12)

## GUT FERMENTATION - HINDGUT ## 
# saving heatmap
ggsave("HINDGUT_HEATMAP.PNG"
       , HINDGUT_HEATMAP
       , height=3, width =23)
# no pca plot to save due to lack of variance between captive vs. wild samples
# saving top 10 upregulated & downregulated pathways
ggsave("HINDGUT_LOG2FC.PNG"
       , HINDGUT_LOG2FC
       , height=6, width =12)
# saving volcano plot
ggsave("HINDGUT_VOLCANO.PNG"
       , HINDGUT_VOLCANO
       , height=6, width =12)

## GUT FERMENTATION - FOREGUT ## 
# saving heatmap
ggsave("FOREGUT_HEATMAP.PNG"
       , FOREGUT_HEATMAP
       , height=13, width =37)
# no pca plot to save due to lack of variance between captive vs. wild samples
# saving top 10 upregulated & downregulated pathways
ggsave("FOREGUT_LOG2FC.PNG"
       , FOREGUT_LOG2FC
       , height=6, width =12)
# saving volcano plot
ggsave("FOREGUT_VOLCANO.PNG"
       , FOREGUT_VOLCANO
       , height=6, width =12)
