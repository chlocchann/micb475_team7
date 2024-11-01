# Meeting Minutes - November 1

## Wendy
* I hate ANCOM :D. It gets confusing with 2+ variables and the ancombc2() function itself takes 20+ minutes to run only to end up with MANY errors. Once you surpass them, there’s the problem of “OTU abundance data must have non-zero dimensions” for a few categories tested, likely meaning that there's no significant groups using this analysis. 
* Phyloseq Q: rarefy_even_depth(sample.size = ?) → data changes depending on this (comparing sample size =1000 vs no restriction)
* DESeq2 (based on sample size =1000) → 
* Tried Jang et al. &  Leon et al.’s code + amend→ both generated similar results, focused on Jang et al.’s 
* Y-axis = captive normalised to wild (so <0 = decrease with adaptive) ; X-axis = phylum
* General conclusion: overall a loss of diversity with captivity (though specific taxonomy groups differ depending on category & group) across mammals
* Note: no significant groups were identified for Omnivores (diet type)
* Note: Bacteroidetes and Firmcutes were present in all remaining comparisons other than Omnivores
* Notice how there are two “positive” groups in Carnivores & Non-fermenter’s Proteobacteria group
* Identified to be the same ASV: "53d42797b8ac1fa460a02c8a2904581e", "f50defe8cf4a491bc5377919cec7e331"
* Samples containing this ASV (all are carnivores + non-fermenters):
* Don't know how to trace back to more data like collector, country, etc :((
### Can we please do DESeq2 and leave ANCOM as future directions :(((

## Tiffany
* Proposal feedback
* In the context of our project, it is important to be specific for…
* Add animal or a synonym before species to clarify that you are not talking about bacterial species.
* Always put “microbial” before “taxa” any time we are referring to microbial taxa as opposed to host taxa
* Intro: adjust first couple paragraphs with more logic, improve flow
* Better build up the need (and benefits/reasons) of investigating differences between captive vs wild animals
### Research objectives
* Further expand why we chose gut fermentation and diet type: Were they not addresses in the discussion or did the authors not continue with analyses exploring these two variables? 
* Research more literature findings to make a hypothesis for the shift in metabolic pathways
* Research aim + rationale
* Work on expanding out aim 2: not laid out as clearly as aim 1 
* Emphasize the link between aim 1 and 2: how does completing aim 1 support analysis towards aim 2? 
* Emphasize how both aims answer our research question 

### Team proposal revisions due Nov 5 → confirm Mya and Christine are leading it
