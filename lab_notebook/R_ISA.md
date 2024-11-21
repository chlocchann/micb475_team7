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
- realized captive/wild wasn't added to gut fermentation ISA and diet types (herbivore / carnivore) wasn't completed
  - clustered GutFermentersHMC and captive_wild variable
  - performed ISA using the clustered variables to achieve the fermentation type along with whether if its captive / wild
  - created tables for each fermentation type + each captive/wild variable
    - combined the tables for fermentation type (eg. hindgut) with the results for both captive and wild variable (eg. made one big table for hindgut captive and wild)
- for the diet types ISA
  - followed the same steps as above with clustering OchHMC and captive_wild, used the same code with a few variable modifications
  - performed ISA using the clustered variables
  - created respective tables for each variable (eg. wild_carnivore, wild_herbivore, etc..) and combined tables for the wild group and the captive group
 - going to double check with TA to check if combining tables is okay
