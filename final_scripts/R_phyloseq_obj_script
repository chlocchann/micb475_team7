#Load packages
library(phyloseq)
library(tidyverse)
library(ape) #importing tree
library(vegan) #for rarecurve fxn


#Load data
meta_data <- read_delim("zoo_metadata.tsv", delim =  "\t")
otu_data <- read_delim("feature-table.txt", delim="\t", skip = 1)
tax_data <- read_delim("taxonomy.tsv", delim = "\t")
phylo_data <- read.tree("tree.nwk")



#### Format OTU table ####

# Replacing the index column with new row names 
otu_mat <- as.matrix(otu_data[,-1]) #remove first column
rownames(otu_mat) <- otu_data$`#OTU ID` #replace w/ OTU IDs

# Creating phyloseq objects
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE) #OTU table
class(OTU)
        # [1] "otu_table"
        # attr(,"package")
        # [1] "phyloseq"


#### Format sample metadata ####

# Save everything except sampleid as new data frame
meta_df <- as.data.frame(meta_data[,-1])

# Make sampleids the rownames
rownames(meta_df)<- meta_data$'#SampleID'

# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(meta_df)
class(SAMP)
        # [1] "sample_data"
        # attr(,"package")
        # [1] "phyloseq"


#### Formatting taxonomy ####

# Convert taxon strings to a table with separate taxa rank columns
tax_mat <- tax_data |>
  select(-Confidence) |> #remove Confidence column
  separate(col=Taxon, sep="; ",
           into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) |>
  as.matrix() # Saving as a matrix

# Save everything except feature IDs 
tax_mat <- tax_mat[,-1]

# Make sampleids the rownames
rownames(tax_mat) <- tax_data$`Feature ID`

# Make taxa table
TAX <- tax_table(tax_mat)
class(TAX)
        # [1] "taxonomyTable"
        # attr(,"package")
        # [1] "phyloseq"


##### Merge all as a phyloseq object #####
zoo_phyloseq <- phyloseq(OTU, SAMP, TAX, phylo_data)

#### Looking at phyloseq object #####
# View components of phyloseq object with the following commands
otu_table(zoo_phyloseq)
sample_data(zoo_phyloseq)
tax_table(zoo_phyloseq)
phy_tree(zoo_phyloseq)



######### ANALYZE ##########
# Remove non-bacterial sequences, if any
zoo_filt <- subset_taxa(zoo_phyloseq,  Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")
# Remove ASVs that have less than 5 counts total
zoo_filt_nolow <- filter_taxa(zoo_filt, function(x) sum(x)>5, prune = TRUE)
# Remove samples with less than 100 reads
zoo_filt_nolow_samps <- prune_samples(sample_sums(zoo_filt_nolow)>100, zoo_filt_nolow)

# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
rarecurve(t(as.data.frame(otu_table(zoo_filt_nolow_samps))), cex=0.1) 
        #transposes OTU column & turn into table
        #creates curve
zoo_rare <- rarefy_even_depth(zoo_filt_nolow_samps, rngseed = 1, sample.size = 10000)
        #rngseed = 1 --> set random # so have reproducibility so sample same way/time
        #1956OTUs were removed because they are no longer  present in any sample after random subsampling
      

##### Saving #####
save(zoo_filt_nolow_samps, file="zoo_filter.RData")
save(zoo_rare, file="zoo_rare_10000.RData")
