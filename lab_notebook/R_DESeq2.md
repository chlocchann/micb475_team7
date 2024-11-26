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

### Oct 29th, 2024 - [WENDY]
- Notice how there are two “positive” groups in Carnivores & Non-fermenter’s Proteobacteria group
- Set up code to extract the ASV identification on R 
- Identified to be the same ASV: "53d42797b8ac1fa460a02c8a2904581e", "f50defe8cf4a491bc5377919cec7e331"
- Set up code to extract the taxonomic information regarding these samples contributing to these ASV
- (Please see meeting minutes or R code for information)

### Oct 30th, 2024 - [WENDY]
- Tried attempting with Ancom, which is more optimal as it has less false positives than DESeq.
- Really struggled with Ancom since we are working with 2+ variables (captivity, diet type/gut fermentation type) and Ancom seems to work mostly with one variable.
- It also took a lot time to run the function (around 20-30 minutes) and the output often returns errors. Will ask TA/group about pursuing with with ANCOM or not.
- DESeq is mostly complete, though concerns may raise about closing 1000 as sample depth form rarefied phyloseq object. Having too many samples also render the data unconclusive (we have 80 000+ samples). Will ask TA/group about sample size.

### Nov 1st, 2024 - [WENDY]
- Decided to continue with DESeq instead and update sample size to 10000 samples. 
- Ended up with many more positive ASVs, though most are still negative.

### Nov 4th, 2023 - [WENDY]
- Decided to blast the previous 2 positive ASVs (from sample size=1000) since there was confusion about sample genus of “Shigella-E. coli”. Confirmed to be E. Coli
- Generated information regarding all the positive ASV (from sample size=10000), and analyzed using excel organizing. Excel also contains taxonomic information of each ASV. Found not many overlaps in ASV between samples. 

### Nov 12th, 2023 - [WENDY]
- The original sample graphs (each phylum composes of different data points representing ASV and bar graphs representing variability in samples under the ASV) may look confusing and disorganized.
- Did not want to do bar graphs because will lose insight into the prevalence of positive/negative valued ASVs (eg. If there are only 1 ASV at value 5, but 100 at value -5, will have a span of -5 to 5 when in reality, most samples are negative)
- Maybe will try violin plots instead, gives the both the benefit of looking at abundance distribution and simplifying representation
- Removed phylums with only one ASV since violin plot cannot plot with just one datapoint (also just one ASV is likely not enough to represent a phylum as significant between captive and adaptive)
- Also tried faceting the graphs together for each variable so can compare phylums as well, especially to see which phylum are commonly significant in the categories of each variable
- Also included code to generate mean and standard deviations of all significant differential ASVs within each phylum in every category/variable


### Nov 26th, 2023 - [WENDY]
- Changed to barplot instead and also recognized it is absurd to have both positive and negative values for a phylum for DESeq plots. Changed up code using template from class, with some adjustments.
- Barplot if remove numerical identification after phylum, collapse graph into having positive and negative spans for some phylum
- Asked Ritu and Evenlyn, decided to keep graph as it is, without excluding identifiers
- Colour coded graphs so anything above log2fold change (captive/wild) of 0 is red and anything below is blue
- Most bars are below 0 (blue), presenting decrease in abundance due to captivity
- Instead decided to generate a code that can summarizes the number of values above and  below 0 for fold change
- Also generated another code that removes the numerical identifier and collapses phylum just by phylum name, and see the number of values above and below 0 for change. This allows us to see which phylum is contributing the most to the overall decrease in phylum abundance (i.e. which phylum contributed to the prevalent observation of bars below 0) and thus, is most affected by captivity.
- Reorganized data above into table view manually for manuscript 


