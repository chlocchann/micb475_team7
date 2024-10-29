### October 25th, 2024 [WENDY]
- Adjusted code in final_scripts/QIIME2_script.sh(renamed table for OTU, adjusted biom convert code, added copy metadata file to export folder code, adjustment of export to local device (added "-r" after scp))) 
- Exported OTU, taxonomy, rooted tree, & metadata files and uploaded to github: qimme2_output_files/exports_for_phyloseq

### Oct 27th, 2024 [WENDY]
- Adapted different DESeq2 to accommodate analyzing different categories of variables (gut fermentation, diet type) between captive and adaptive:
- Course code: results in a massive complement of groups on X against log2change that is hard to distinguish letters by eye or conclude due to massive amount of positive data; tried different levels: genus, family, order, and phylum but illegibility and inconclusiveness still exist
- Jang et al. (Sex Influences Gut Microbial Composition in Mice with Familial Dysautonomia but is not the Primary Determinant of Microbial Functional Diversity, UJEMI 09 2024): manuscript shows 10 phylum groups across X axis; adapted code to our scenario but still resulted in same situation as above
- Leon et al. (Captivity plays a role in shaping the gut microbiome of social mammals, UJEMI, 09 2022): manuscript shows 4 phylum groups across X axis; adapted code to our scenario but showed both positive and negative values across a phylum group, which did not appear in previous examples nor does it make sense
- Possible error in massive amount of data due to used filtered phyloseq obj data instead of rarefied (took too long to run on R so decided to go test filtered first)

### Oct 28th, 2024 [WENDY]
- Indeed error was due to unrarefied data
- Tried with rarefied data (depth =1000) and results were much better:
- Leon et al. shows 2, 5, 6 phylum groups respective for foregut, hindgut, non-fermentation. First two shows all negative log2change (x-axis) values while last one (non-fermentor) has the span positive and negative values again (Proteobacteria).
- Jang et al. looks more realistic in terms of span but found the x groups are not well classified (eg. P__Bacteroidota and P__Bacteroidota.1 are considered separate groups even through they are the same)
- Rearranged Jang et al. (added "mutate(Phylum_Clean = gsub("p__|\\..*", "", Phylum)) %>%” so can better group), result also 5 groups like Leon et. al by log fold change is /10 in Jang et al.’s case
- After resolving this problem, tried looking at non-fermentor and found Proteobacteria to have quite a large span by individual points, likely reason why span ositive and negative values in Leon et al.. Also note, non-fermentor only had 5 groups instead of 6.

