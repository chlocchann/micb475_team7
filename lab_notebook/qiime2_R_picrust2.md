### November 4, 2024 [CHLOE]
* Added PICRUSt2 pipeline to existing R script
* Ran PICRUSt2 pipeline and transferred files to local device
* Edited "pathway_abundance.tsv" such that first row contained #OTU IDs
* Uploaded to RStudio alongside metadata file
* Adapted code from "picrust_analysis.R" taught in lecture (up until # generate pathway PCA plot)
* Ran into issues with correlating captivity status (captive_wild column) with diet type (OchHMC) and gut fermentation type (GutFermentersHMC)
  * Current code takes in one metadata column where 2 variables are compared -> must compile captivity status and metadata columns of interest into one column? 

### November 5, 2024 [CHLOE]
* Appended " - CAPTIVE" and " - WILD" to diet type (OchHMC) and gut fermentation type (GutFermentersHMC)
* Tested code from "picrust_analysis.R" taught in lecture (up until # generate pathway PCA plot)
  * Edited pipeline such that only one variable comparing captive vs wild status was being compared for each comparison group
  * Only 8 samples remained for omnivore group after filtering -> 6 captive, 2 wild
    * ERROR: heatmap only shows captive omnivores despite indicating that wild samples are shown ?
  * ERROR: PCA pathway plots for omnivore & herbivore: Error in prcomp.default(t(abundance), center = TRUE, scale = TRUE :cannot rescale a constant/zero column to unit variance
    * Troubleshooting suggests that there is no variance between captive and wild samples for both of the above
  * ERROR: PCA pathway plots for hindgut and foregut fermenters: Error in prcomp.default(t(abundance), center = TRUE, scale = TRUE :cannot rescale a constant/zero column to unit variance
    * Troubleshooting suggests that there is no variance between captive and wild samples for both of the above