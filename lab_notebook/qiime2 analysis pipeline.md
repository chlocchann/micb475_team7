### October 10, 2024 [CHLOE]
- Created working directory under /data -> named team7_captivevswild
- Copied seqs, manifest & metadata from /mnt/datasets/project_2/zoo into working directory (NOTE: both manifest & metadata are .txt files)
- Created zoo_demux detached screen for demultiplexing & ran overnight

### October 11, 2024 [CHLOE]
- Created demux.qzv visualization file & transferred to local device
- Visualized with qiime2 visualization tool & decided on truncation length of 233bp
- Created zoo_DADA2 detached screen for determining ASVs & ran overnight

### October 12, 2024 [CHLOE]
- Created denoise_stats.qzv, zoo_table.qzv, and rep_seqs.qzv visualization files & transferred to local device
- Created zoo_silva detached screen for training classifier & ran overnight
- Extracted V4 16s rRNA regions with the 515F (Parada)-806R (Appril) V4 amplification primers using Silva Database as reference

### October 13, 2024 [CHLOE]
- Trained classifier with trimmed ref-seq file in zoo_silva detached screen and ran overnight
  
### October 14, 2024 [CHLOE]
- Assigned taxonomy to reads using trained classifier
- Created taxonomy.qzv visualization file and transferrred to local device
- Filtered out mitochondrial & chloroplast samples
- Created zoo_filtered_table.qzv visualization file and transferred to local device
- Made taxonomy bar plots
- Created taxa_bar_plots.qzv visualization file and transferred to local device
- Created zoo_phylo detached screen for phylogenetic diversity analyses
- Generated tree for phylogenetic diversity analyses and ran overnight
  
### October 15, 2024 [CHLOE]
- In zoo_phylo detached screen, ran alpha-rarefaction from morning till evening
- Error message detected: Plugin error from diversity: [Errno 2] No such file or directory: '/home/qiime2/TMP/qiime2-temp-grwe7aai6/shannon-21-1_PopulationDensity_n/km2.jsonp' Debug into has been saved to /home/qiime2/TMP/qiime2-q2cli-err-zkpzidln.log
- Emailed Ritu for assistance -> recommended converting metadata.txt file to a metadata.tsv file to see if remedies issue
- Ran alpha-rarefaction with metadata.tsv file overnight
  
### October 16, 2024 [CHLOE]
- Error message persisted (same as above)
- Emailed Ritu to update -> said mentions Shannon's diversity metric
- Noticed that it referenced particular metadata column "21-1_PopulationDensity_n/km2" -> hypothesized may have something to do with particular column and opted to remove and retry
- Ran entire qiime2 workflow up to training clasifier throughout the day (needed to extract 16s rRNA overnight)
  
### October 17, 2024 [CHLOE]
- Ran qiime2 workflow from training classifier to alpha-rarefaction throughout the day
- Ran alpha-rarefaction with edited metadata.tsv file (excludes 21-1_PopulationDensity_n/km2 column) and ran overnight

### October 18, 2024 [CHLOE]
- Error message updated: Plugin error from diversity: [Errno 2] No such file or directory: '/home/qiime2/TMP/qiime2-temp-grwe7aai6/shannon-HuPopDen_Min_n/km2.jsonp' Debug into has been saved to /home/qiime2/TMP/qiime2-q2cli-err-jhju_08a.log
- Browsed qiime2 forums to troubleshoot -> found similar [issue](https://forum.qiime2.org/t/plugin-error-from-diversity-errno-2-no-such-file-or-directory/18892) indicating problem lies in presence of "/" in column names -> removed in excel and exported as .tsv
- Reuploaded new zoo_metadata.tsv file and re-ran alpha-rarefaction -> SUCCESS!!
