# October 18 Meeting Minutes

## Tiff’s questions
* Variables used: we discussed last time to do genus, diet type, and gut fermentation but not sure about genus anymore → too many categories to compare, is it a useful variable?
* Do we have to make a hypothesis in Research objective section? What can it be?
* Hypothesize that there will be differences? 
* should we not do within categories? Comparisons→
1) wild vs captive non-fermenters
2) wild vs captive hind-fermenters
3) wild vs captive fore-fermenters
4) within wild fermenters: fore vs hind
5) within captive fermenters: fore vs hind
Should we group them together? So it’s 4 categories



## Wendy’s questions
* Where do we draw the line between background and research objective?
* It sounds repetitive in some parts: 
* Ie. do we introduce “our study aims to do this__” in the intro? 

## QIIME2 Pipeline Update
* Had issues with running alpha-rarefaction step of the pipeline -> pinpointed to presence of “/” in certain metadata column names -> removed and worked!
* Does removal of “/” from metadata column names need to be “codified”? (ie. does it need to be done in R Studio using something like rename() or can it just be removed in Excel and a comment made a the top of the QIIME2 pipeline)
