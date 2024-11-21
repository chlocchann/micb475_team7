### NOV 13, 2024 [Christine]
- loaded rarefied phyloseq object [zoo_rare]
- used code learned in class: aggregated taxa to the species level by using tax_glom then transformed to relative abundance using transform_sample_counts
- ran the indicator species analysis, adapted the code from class [multipatt] and adjusted the variables to use the relative abundance as well as 'GutFermentersHMC'

### NOV 14, 2024 [Christine]
- extracted taxonomic table from zoo_rare into a new dataframe and added an ASV column for filtering purposes later on
- merged the taxonomy table with the phyloseq object and filtered out non-significant values (p-value>0.05) and indicator values (stat<0.80)
- for each fermentation type (foregut, hindgut, non-fermenters), created a table that includes the ASV, Phylum, p-value and indicator value (stat)

### NOV 20, 2024 [Christine]
- tried creating a table that includes all of the taxonomic levels, but looked overly confusing (species column included NAs which made it more confusing)
- just kept the original table for each fermentation type 

