# Meeting Minutes - October 11

# Determine our research question

# Project Proposal Outline
* Discuss tasks/roles 
* Discuss timeline

# Data Overview/Data Wrangling Section
* CLARIFICATION: Do we not need to include information on how we’re planning to process our data using QIIME2 in the proposal if the Data Wrangling rubric replaces the Data Overview one?
* Data Wrangling focuses on defining metadata to better suit purposes of modelling -> ex. Changing numericals to categoricals
* Data Overview focuses on entire QIIME2 pipeline -> according to previous modelling projects, still need to go through QIIME2; should we include both the data overview and data wrangling information in our proposal?

## DENOISING
![Screenshot 2024-10-10 at 1 50 43 PM](https://github.com/user-attachments/assets/5c8e74c1-b63b-48e9-9911-867e41dad5fc)
* Thinking of trimming at 233 bp -> includes low quality read at 206bp, but average of reads is still > Phred Score = 30; worth it to preserve base information? 
* Original paper trimmed at 100bp
* Unsure what lack of boxes means on Interactive Quality Plot? (sprinkled in starting at 183bp)

## TRAINING A CLASSIFIER
* Original paper assigned taxonomy using the RDP classifier and the Greengenes August 2013 release as the reference database.
* DNA amplification performed on V4 region of 16s rRNA gene -> 515F (Caporaso)–806R (Caporaso) primers -> FWD:GTGCCAGCMGCCGCGGTAA;REV:GGACTACHVGGGTWTCTAAT
* Updated V4 primers (forward, 2016; reverse, 2015): 515F (Parada)–806R (Apprill)
FWD:GTGYCAGCMGCCGCGGTAA; REV:GGACTACNVGGGTWTCTAAT

## DIVERSITY ANALYSES
* Original paper used Shannons (Alpha) & Bray-Curtis (Beta)

# Tiff’s Thoughts 
## Determine which variables to use/filter out 
- The paper already studied which host factors contribute to difference in gut microbial diversity in wild vs captive animals (diet type, fermentation type, genus)
- If we’re still answering the question: what metadata characteristics are associated with increased/decreased microbial diversity in the mammalian gut microbiome? →  which variables should we use?
### Metadata categories
- Animal type
- Location (ie. from zoo/park)
- The geographic location of the samples have different environmental factors that can affect gut microbiome diversity.
- Country
- Differences in management practices in captivity, which could affect gut microbiome diversity
- Taxonomy (ie. phylum, class, order, family, genus, species) → they already did genus so what reason do we have to do family? Species is too specific..
- Activity (ie. crepuscular, diurnal, nocturnal)
- Activity levels can impact dietary habits therefore, gut microbiome composition. Is there a behavioural pattern here that correlates with diversity level?
- Diet → they already did it
- Trophic level
- Conservation status → they showed not significant 

## Go through our aims, do we need to hit 4?

## Research Objective VS aim/rationale section
Is the first section background information about why for ex, we need to research which variables are the most likely to be associated. Then the next section is how we plan to identify those variables?
