# Meeting Minutes - November 22

## Chloe (PICRUSt2)
* All graphs/plots have been generated -> PICRUSt2 Graphs & Plots - MICB 475 (in Google Drive) 
	* IMO, top upregulated/downregulated bar plots & volcano plots should go in results, while heatmaps & pca should go in supplementary
* Entire PICRUSt2 code (incl edited DESEq2 function) is in github alongside updated lab notebook
* Results section is complete -> unable to make inferences about particular pathways, so just pointed out general trends (which comparisons had more/less pathway abundance changes)
* Wrote code for metadata wrangling (column name changes & appending captivity status)
	* Should methods section state column name changes? No right?
* PICRUSt2 methods -> do we cite TA who wrote the custom DESEq2 function?
* Shown on the bar graphs is what’s downregulated (negative values) for the captive group. Carnivore and omnivore is all downregulated. Data matches with ISA data → loss of metabolic function 
* Put the bar graphs as supplementary, use the volcano plots for results 

## Wendy (DESEq2)
* Updates from Evelyn on if did organising on excel is ok vs req design code for that?
* Deseq violin plots MICB475_DESeq2 PNGs.pptx: removed ones with just one ASV since violin plot cannot plot with just datapoint
	* For explanation: both cannot plot + just one ASV not enough to represent phylum
	* “Phyla with only one significant differential ASV were excluded to improve plot clarity and representability of abundance changes at the phylum level”
	* Also tried facetting so will able to see which ones are avail in all categories
	* Mean and variance summarised
* Positive ASVs also summarised
* Discussion (Deseq)
	* (results) Decreased ASVs abundance in captive animals compared to wild animal at the phylum level in all categories within A) diet type and B) gut fermentation, 
	* (discussion) meaning decreased in microbial diversity, as of expectations compared to McKenzie et al.’s (beta/alpha diversity analysis). This was also observed in (previous literature) …
	* Addressing the lack of resolution into which microbial communities were significantly decreased, our study sorted differential ASVs in wild versus captive of each category within the variables at the phylum level.
	* Specifically, we found that:
		* Bacteriodota and Firmcutes were present and majorly decreased in captive (vs wild)  in all graphs → look into its specific functions to explain why it decreases
		* Verrucomicrobiota present in all gut fermentation type and decreased
		* Optional: talk about how Proteobacteria and Spirochaetota have positive mean ASV in many situations
		* Herbivore, Proteobacteria: 0.23503389 (mean), 5.865863752 (sd)
		* Carnivore, Proteobacteria: 261349899 (mean), 3.391782587 (sd)
		* Foregut, Spirochaetota: 1.657513507 (mean), 5.172504123 (sd)
		* Non-fermentors, Proteobacteria: 0.032465117 (mean), 4.778983847 (sd)
		* Non-fermentors, Spirochaetota: 2.717835307 (mean), 2.976781032 (sd)
		* Also please mention something about how although mean is positive, standard deviation is also much higher than other phylums in category to help support major result
	* Conclude how these seem to be generally key players affected by captivity (Bacteriodota + Firmcutes for all and Verrucomicrobiota for gut fermentation type) across all mammals
	* Mya/Tiffany: please make this sound better + look into highlighted portions (observed in (previous literature) & Bacteriodota & Firmcutes & Verrucomicrobiota & maybe Proteobacteria and Spirochaetota) + connect w/ Chloe’s part for their potential roles in the identified fxnal pathways
* Limitations (Deseq)
	* Only performed on rarefied samples (sample size of 10000) when we have 83,290 samples. Did so for better visualisation because of the large number of samples. Expect probably will not change major trends if increase but would increase the number of ASVs and provide more insight into less represented phylums, such as the ones removed due to having only one significantly differential ASV between captive and wild.
	* Not all significantly differential ASVs were decreased in captive compared to wild
* Future directions (Deseq)
	* Significantly differential ASVs increased (log2fold >0) in captivity compared to wild cha identified and calculated number of times it showed up across the categories in the variables. Data is in supplementary, can further perform more analysis on them to identify trend such as specific animal, (etc categories in the metadata excel)
	* Looking into specific categories within each variable because right now we are just summarising
	* While DESeq provides good insight, ANCOM is also a possible alternative as it captures less false positives 

## Christine (ISA)
* Have about 5 different tables for each respective fermentation type (hindgut, foregut, non-fermenters) and diet types (carnivore, herbivore)
* The generated tables → ISA tables , ( NOT FINAL; gotta double check with the questions below); had a hard time putting the table together, not sure how to ‘nicely’ show it
* The ISA gut fermentation results → ISA GUT results  &  ISA diet type results → ISA DIET results
	** the cells are coloured for my own reference, will revert to original if needed later on
* QUESTIONS: 
	* Should I include all of the taxonomic levels into each table, OR should I just include the ASVs and Phylum?
* Formatting tables
	* Should I remove NAs from species? I didn’t remove it initially since high IV were on species that were ‘NA’
*TO-DO:
	* Paste ISA results portion after meeting (need clarification from TA on a couple things)

* We want to know how many are unique in each category, not the one that are shared between categories → should do core microbiome analysis
	* Perhaps present this as the table of ISA not the exhaustive list
	* See a general pattern, wild has more indicator species in all the categories 
	* Available on github, don’t put it on supplementary 
	* Wild and captive at the top of table 

## Tiffany
* Abstract, intro, and methods draft is complete, waiting for remaining parts of manuscript before finishing up 
* QUESTIONS:
	* MICB_475 DESeq2 PNGs.pptx (linked in Google Drive) 
	* p.2 top violin plots, is one of them supposed to be carnivores? 
	* p.2 violin plots, are we using the top ones or bottom? From my understanding, the bottom ones are to show the same phylums left-right for easy comparison, but would it be confusing if some of the phylums are empty lol
	* some of the lines are straight/do not show density (ie. non-fermenters): how do we read that? Is it accurate to say this in results:
	* Overall, our results showed that most phyla were reduced in captive animals compared to wild animals, as evidenced by negative fold change values relative to the wild group, regardless of category or variable. 
* PICRUSt2 Graphs & Plots - MICB 475
	* Bar graphs: would it be possible for the gradient displaying the significance to be consistent across different plots for that variable?
	* Bar graphs: should we have it in the same scale across categories?
	* Results: I think we should be more specific? Instead of saying  all diet types had functionally distinct metabolic pathway changes in captive animals compared to their wild counterparts → what does the trend look like? Which pathways or what majority were up vs down regulated? Is the trend consistent across categories or are they different, how are they different?

## Action items
* do core microbiome analysis
* do core microbiome to tell you how much is unique and how much is common 
* ISA assign "association" with different groups (assess prevalence and abundance)
* do the 6 venn diagrams but put it as a table
* can combine it into 1 table for ISA + core microbiome 
* Fix up DESEq: should only have 1 fold change for each phylum 
* Discussion: pick the top in each analysis, talk about it a bit but mostly focus on the larger trends 
