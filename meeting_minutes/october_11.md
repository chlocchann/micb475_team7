# Meeting Minutes - October 11

# INTRODUCTION SECTION
What is the knowledge gap?
The paper’s rationale was that few studies have yet to address how changes in living environments can influence mammalian gut microbiome, and for the current studies, they target few species → so they studied effect of captivity on the mammalian gut microbiome → Is our rationale just expanding on that, and finding more associated variables?
They didn’t do a functional analysis 
We are looking at functional analysis + Which ones are enriched in the significant variables + what is the distribution of these microbes 
Have some evidence that variables XYZ have taxonomic differences, analyse the functional pathways that are enriched 
But we still need to repeat taxonomic from them as well 
So are we still doing modelling or focusing on the variables that have already shown significance? 
Leaning towards focusing on the significant variables: host taxonomy (genus), diet type, and gut fermentation type → can skip finding variables (step 1)
Taxonomic analysis: DEseq2
Look at abundance 
Wild vs captive
Within wild
Within captive
Functional analysis
How does X, Y Z variable influence metabolic pathways in wild vs captive mammals?
Within wild animals, does X, Y Z variable influence metabolic pathways?
Within captive animals, does X, Y Z variable influence metabolic pathways?

# RESEARCH OBJECTIVES SECTION
## confirm research question: can it be broad like… what variables are associated with changes in microbial diversity in the mammalian gut microbiome of wild vs captive mammals? 
Is this sufficient to explain why we’re including the variables already tested: checking if our pipeline produces the same results as the paper for host variables used

## Do we go need to go through each of the metadata variables (table below)  selected and explain why we are including them?
hypotheses: would it be explaining “we hypothesize that variable X will be associated with changes in gut microbial diversity in captive vs wild animals” or would it be like “we expect variable X to be the most significant variable associated…”

## Decide what metadata variables to use zoo_metadata.xlsx
Should we use numerical variables (ie. body  mass) and change it into categorical?
Metadata categories + rationale for keeping them (questions in red)
Category
Rationale
Location (ie. from zoo/park) 
The geographic location of the samples have different environmental factors that can affect gut microbiome diversity. 
Country 
Differences in management practices in captivity, which could affect gut microbiome diversity
Lifestage

Animal type (ie. Common name) 

Taxonomy (KPCOFGS)
Note: paper already did genus so what reason do we have to do family or higher orders?

Would species be too specific?
Activity (ie. crepuscular, diurnal, nocturnal)
Activity levels can impact dietary habits therefore, gut microbiome composition. Basically answering the question, is there a behavioural pattern here that correlates with diversity level?
Body mass 
Numerical so we need to convert it to category → how do we split and is it worth splitting?
Diet
Note: paper did diet type by herbivore, carnivore, omnivore, we can use all the Eltonian trait diet categories?
Diet Breadth
Would this be worth doing?  Since it is numerical, would each number be a category (1-5)?
Habitat Breadth
Category #1-3? Can answer whether higher habitat breadth is correlated to changes in microbial diversity 
Max longevity?
Social group size 
Sociality 
Trophic level 
Conservation Status
They showed not significant but we can reaffirm
Gut fermentation 
They showed significant but we can reaffirm

# EXPERIMENTAL AIMS SECTION
## Are we missing anything else?
Aim 1: Identify which variables are significantly associated with gut microbiome diversity differences in wild vs captive mammalians.
Bray-curtis community dissimilarity between wild + captive
Rationale: provide insight into which variables are associated with significant changes to the gut microbial community composition between wild and captive animals, which will allow us to further investigate the gut microbial dynamics in each variable (aim 2) and metabolic pathways affected (aim 3)

Aim 2: Examine microbial taxa and alpha and beta diversity on variable(s) of interest that significantly correlate to differences in gut microbiome diversity in wild vs captive mammalians. 
For example, in the types of gut fermentation, are there differences in microbial diversity + microbial taxa?
What is the difference between alpha and beta university in this context?
If the variable is gut fermentation type, and we’re comparing wild vs captive mammalians, would alpha be looking at variations of microbes within each of hindgut, midgut? And beta looks at how different are the microbial communities from the hindgut vs midgut?
In the variable of interest, use DESeq to assess differences in microbial composition at the ASV level, then we can look at those ASVs (what taxonomic group do they belong to?)
Rationale: help us determine how gut microbiome diversity varies in each variable of interest? Help us determine if certain microbial taxa are significantly correlated to differences in gut microbiome diversity in wild vs captive mammals. To investigate the microbial structures within the variable of interest

Aim 3: Functional analysis of variable(s) of interest. 
Are we predicting metabolic pathways based on significant microbial taxa from aim 2? Or do we just look into the signfiicant variables found in aim 1? 
predict what metabolic pathways are enriched or depleted in relation to the significant variables linked to gut microbiome diversity differences in wild vs captive mammals
Rationale: 

# Data Overview/Data Wrangling Section
## CLARIFICATION: Do we not need to include information on how we’re planning to process our data using QIIME2 in the proposal if the Data Wrangling rubric replaces the Data Overview one? 
* Data Wrangling focuses on defining metadata to better suit purposes of modelling -> ex. Changing numericals to categoricals
* Data Overview focuses on entire QIIME2 pipeline -> according to previous modelling projects, still need to go through QIIME2; should we include both the data overview and data wrangling information in our proposal?

## DENOISING

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


# POSSIBLE APPROACH
## ”Sex Influences Gut Microbial Composition in Mice with Familial Dysautonomia but is not the Primary Determinant of Microbial Functional Diversity”:
QIIME2 processing, DAD2 denoising and clustering, taxonomic classification of the ASV w/ SILVA (mod 6), filtered out chloroplast & mitochondria
R tidyverse phyloseq object on result of (1), removed low ASV reads = unrarefied phyloseq file
Created rarefied phyloseq file (changed rarefied sampling depth) on result of (2) = rarefied phyloseq file
Filtered & retained variables on (3, rarefied), Bray-Curtis distance (basis for variation calculation for PERMANOVA), PERMANOVA (multivariate analysis, P-value, R^2), identify significant variables → OUR AIM1: VARIABLE SIGNIF
		
Create DESeq file on result of (2, unrarefied) via DESeq2, normalised group 2 to group 1, volcano plot (P-value), identify significant and non-significant ASVs → OUR AIM2: TAXONOMIC ANALYSIS
	
Also did smt called Indicator Species Analysis (ISA) for taxonomic analysis:

PICRUSt2 on QIIME2 (CANVAS mod) on removed ASV <5 version of (1), mapped predicted gene family abundances to known metabolic pathways into pathway abundance table, table  converted to a human-readable file compare metabolic pathways, abundance table filtered to contain samples that existed in metadata of the rarified phyloseq object (3), pathway identifiers annotated w/ MetaCyc pathway = relative abundance of significant metabolic pathway file
Create Principal Component Analysis (PCA) plot, identify diverse metabolic pathway abundance across samples btw group1 & group2 → OUR AIM3: FXNAL ANALYSIS
(Fig 4A, identified distribution, and commented on overlap (no clear separation in the abundance of metabolic pathway btw groups)
Kinda confusing, doesn't really show too much data tbh

Conducted heat wave using same annotated file of (6), identify upregulated pathways in group 1 compared to group 2 → OUR AIM3: FXNAL ANALYSIS
(Fig 4B)
Spearman’s rank correlation test on scatter plot created from calculated relative abundance (%) at the phylum level using (2, unrarefied) w/ dplyr (see 10a.) + result of (6)

## Other things:
Taxonomic composition analysis on (2, unrarefied) w/ calculating phylum level using dplyr and ggplot histograms → Fig 2
Alpha and beta analysis → Supplementary
