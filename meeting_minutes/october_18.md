# Meeting Minutes - October 18

## Tiff’s questions
* Variables used: we discussed last time to do genus, diet type, and gut fermentation but not sure about genus anymore → too many categories to compare, is it a useful variable?
  * Can leave out genus, focus on other variables (diet type, gut fermentation)
  * Do we have to make a hypothesis in Research objective section? What can it be?
* Hypothesize that there will be taxonomic differences and functional differences: based on? 
* Hypothesize that functional pathways associated with specific bacterial taxa will be the ones that are enriched
  * Can combine taxonomic + functional differences into 1, don’t need to be too specific for objective portion
* What are the different levels of comparisons 
  * Aim 1: Taxonomic comparison ISA and DESEq2
    * 1A first variable 
    * 1B second variable 
  * Aim 2: Functional comparison PICRUrst
    * 2A first variable 
    * 2B second variable 
    * Comparisons→
      * 1) wild vs captive non-fermenters
      * 2) wild vs captive hind-fermenters
      * 3) wild vs captive fore-fermenters
## Wendy’s questions
* Where do we draw the line between background and research objective?
  * It sounds repetitive in some parts: 
    * Ie. do we introduce “our study aims to do this__” in the intro? 

## QIIME2 Pipeline Update (Chloe)
* Had issues with running alpha-rarefaction step of the pipeline -> pinpointed to presence of “/” in certain metadata column names -> removed and worked!
* Does removal of “/” from metadata column names need to be “codified”? (ie. does it need to be done in R Studio using something like rename() or can it just be removed in Excel and a comment made a the top of the QIIME2 pipeline) -> WORRY ABOUT IT LATER, SHOULD BE FINE FOR NOW
* Reference to QIIM2, R packages, SILVA, primers (Earth smt org)
* Decide sampling depth based on fermentation type and diet
