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
* Adapted code from "picrust_analysis.R" taught in lecture (full code)
  * Edited pipeline such that only one variable comparing captive vs wild status was being compared for each comparison group
