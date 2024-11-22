# Meeting Minutes - November 22

## Chloe (PICRUSt2)
* All graphs/plots have been generated -> PICRUSt2 Graphs & Plots - MICB 475
	** Adjusted code based on presence of upregulated/downregulated pathways
		*** If only downregulated, only filtered for top 10 downregulated
		*** If less than 10 up/downregulated, did not filter further
** IMO, top upregulated/downregulated bar plots & volcano plots should go in results, while heatmaps & pca should go in supplementary
		*** Entire PICRUSt2 code (incl edited DESEq2 function) is in github alongside updated lab notebook
* Results section is complete -> unable to make inferences about particular pathways, so just pointed out general trends (which comparisons had more/less pathway abundance changes)
* Wrote code for metadata wrangling (column name changes & appending captivity status)
	* Should methods section state column name changes? No right?
* PICRUSt2 methods -> do we cite TA who wrote the custom DESEq2 function?

## Wendy
* updates from Evelyn on if did organising on excel is ok vs req design code for that?
* Deseq violin plots MICB_475 DESeq2 PNGs.pptx: removed ones with just one ASV since violin plot cannot plot with just datapoint
For explanation: both cannot plot + just one ASV not enough to represent phylum
“Phyla with only one significant differential ASV were excluded to improve plot clarity and representability of abundance changes at the phylum level”
Also tried facetting so will able to see which ones are avail in all categories
* discussion (Deseq)
(results) Decreased ASVs abundance in captive animals compared to wild animal at the phylum level in all categories within A) diet type and B) gut fermentation, 
(discussion) meaning decreased in microbial diversity, as of expectations compared to McKenzie et al.’s (beta/alpha diversity analysis). This was also observed in (previous literature) …
Addressing the lack of resolution into which microbial communities were significantly decreased, our study sorted differential ASVs in wild versus captive of each category within the variables at the phylum level.
Specifically, we found that:
Bacteriodota and Firmcutes were present and majorly decreased in captive (vs wild)  in all graphs → look into its specific functions to explain why it decreases
Verrucomicrobiota present in all gut fermentation type and decreased
Optional: talk about how Proteobacteria and Spirochaetota have positive mean ASV in many situations
Herbivore, Proteobacteria: 0.23503389 (mean), 5.865863752 (sd)
Carnivore, Proteobacteria: 261349899 (mean), 3.391782587 (sd)
Foregut, Spirochaetota: 1.657513507 (mean), 5.172504123 (sd)
Non-fermentors, Proteobacteria: 0.032465117 (mean), 4.778983847 (sd)
Non-fermentors, Spirochaetota: 2.717835307 (mean), 2.976781032 (sd)
Also please mention something about how although mean is positive, standard deviation is also much higher than other phylums in category to help support major result
Conclude how these seem to be generally key players affected by captivity (Bacteriodota + Firmcutes for all and Verrucomicrobiota for gut fermentation type) across all mammals
Mya/Tiffany: please make this sound better + look into highlighted portions (observed in (previous literature) & Bacteriodota & Firmcutes & Verrucomicrobiota & maybe Proteobacteria and Spirochaetota) + connect w/ Chloe’s part for their potential roles in the identified fxnal pathways
* limitations (Deseq)
Only performed on rarefied samples (sample size of 10000) when we have 83,290 samples. Did so for better visualisation because of the large number of samples. Expect probably will not change major trends if increase but would increase the number of ASVs and provide more insight into less represented phylums, such as the ones removed due to having only one significantly differential ASV between captive and wild.
Not all significantly differential ASVs were decreased in captive compared to wild
* future directions (Deseq)
Significantly differential ASVs increased (log2fold >0) in captivity compared to wild cha identified and calculated number of times it showed up across the categories in the variables. Data is in supplementary, can further perform more analysis on them to identify trend such as specific animal, (etc categories in the metadata excel)
Looking into specific categories within each variable because right now we are just summarising
While DESeq provides good insight, ANCOM is also a possible alternative as it captures less false positives 


## Christine (ISA)
* have about 5 different tables for each respective fermentation type (hindgut, foregut, non-fermenters) and diet types (carnivore, herbivore)
* the generated tables → ISA tables , ( NOT FINAL; gotta double check with the questions below)
* QUESTIONS: 
** Should I include all of the taxonomic levels into each table, OR should I just include the ASVs and Phylum?
** combining tables?
** if clustering variables (eg. gut ferm + captivity, and OchHMC + captivity) is okay??
** I was thinking about clustering all 3 variables together, but not sure if that would be very useful

##Tiffany
* abstract, intro, and methods draft is complete, waiting for remaining parts of manuscript before finishing up 
* powerpoint draft
