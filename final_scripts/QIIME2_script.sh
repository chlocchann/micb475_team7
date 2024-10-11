# ! /bin/bash

## LOGIN CREDENTIALS ##

ssh root@10.19.139.186
# password: Biome1594

## CREATE A WORKING DIRECTORY FOR TEAM 7 ##

# making the working directory
cd /data
mkdir team7_captivevswild
cd team7_captivevswild

# copy zoo data into working directory
cp -r /mnt/datasets/project_2/zoo/seqs .
cp /mnt/datasets/project_2/zoo/zoo_manifest.txt .
cp /mnt/datasets/project_2/zoo/zoo_metadata.txt .

## IMPORTING AND DEMULTIPLEXING USING MANIFEST ## 

# create a detached screen session & activate QIIME2 to run overnight
screen -S zoo_demux
conda activate qiime2-2023.7

# importing & demultiplexing
 qiime tools import \
--type "SampleData[SequencesWithQuality]" \
--input-format SingleEndFastqManifestPhred33V2 \
--input-path ./zoo_manifest.txt \
--output-path ./demux_seqs.qza

## VISUALIZING DEMULTIPLEXED FILES ## 

# creating .qzv
qiime demux summarize \
  --i-data demux_seqs.qza \
  --o-visualization demux.qzv

# transferring demux file to local device
scp root@10.19.139.186:/data/team7_captivevswild/demux.qzv .

## DETERMINING & VISUALIZING ASVs ##

# determining ASVs with DADA2
# based on qiime2 visualization of demux.qza, truncation length of 233bp was selected
qiime dada2 denoise-single \
  --i-demultiplexed-seqs demux_seqs.qza \
  --p-trim-left 0 \
  --p-trunc-len 233 \
  --o-representative-sequences rep_seqs.qza \
  --o-table zoo_table.qza \
  --o-denoising-stats denoise_stats.qza

# visualizing DADA2 stats
qiime metadata tabulate \
  --m-input-file denoise_stats.qza \
  --o-visualization denoise_stats.qzv

qiime feature-table summarize \
  --i-table zoo_table.qza \
  --o-visualization zoo_table.qzv \
  --m-sample-metadata-file /data/team7_captivevswild/zoo_metadata.txt
  
qiime feature-table tabulate-seqs \
  --i-data rep_seqs.qza \
  --o-visualization rep_seqs.qzv

# transferring visualization files to local device
scp root@10.19.139.186:/data/team7_captivevswild/denoise_stats.qzv .
scp root@10.19.139.186:/data/team7_captivevswild/zoo_table.qzv .
scp root@10.19.139.186:/data/team7_captivevswild/rep_seqs.qzv .

## TRAINING & USING CLASSIFIER FOR TAXONOMIC ANALYSIS ##
# extracting V4 16s rRNA region from Silva Database
# used 515F (Parada)–806R (Apprill) V4 amplification primers
qiime feature-classifier extract-reads \
  --i-sequences /mnt/datasets/silva_ref_files/silva-138-99-seqs.qza \
  --p-f-primer GTGYCAGCMGCCGCGGTAA \
  --p-r-primer GGACTACNVGGGTWTCTAAT \
  --p-trunc-len 250 \
  --o-reads ref-seqs-trimmed.qza

# training classifier with trimmed ref-seq file
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ref-seqs-trimmed.qza \
  --i-reference-taxonomy /mnt/datasets/silva_ref_files/silva-138-99-tax.qza \
  --o-classifier classifier.qza

# assigning taxonomy to reads using trained classifier
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

# visualizing taxonomy file
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

## FILTERING & TAXONOMIC ANALYSIS ##

# filtering out mitochondrial and chloroplast samples
qiime taxa filter-table \
  --i-table zoo_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table zoo_filtered_table.qza

qiime feature-table summarize \
  --i-table zoo_filtered_table.qza \
  --o-visualization zoo_filtered_table.qzv \
  --m-sample-metadata-file /data/team7_captivevswild/zoo_metadata.txt

# taxonomy barplots
qiime taxa barplot \
  --i-table zoo_filtered_table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file /data/team7_captivevswild/zoo_metadata.txt \
  --o-visualization taxa_bar_plots.qzv

# transferring taxa bar plot visualization file to local device
scp root@10.19.139.186:/data/team7_captivevswild/taxa_bar_plots.qzv .

## PHYLOGENETIC DIVERSITY ANALYSES ##
# [INSERT]

## EXPORTING OTU, TAXONOMY & ROOTED TREE FILES ##
# creating directory for export
mkdir team7_exports

# OTU table export
qiime tools export \
  --input-path /data/team7_captivevswild/table.qza \
  --output-path /data/team7_captivevswild/team7_exports

biom convert -i feature-table.biom --to-tsv -o feature-table.txt

# taxonomy export
qiime tools export \
  --input-path /data/team7_captivevswild/taxonomy.qza \
  --output-path /data/team7_captivevswild/team7_exports

# rooted tree export
qiime tools export \
  --input-path /data/team7_captivevswild/rooted-tree.qza \
  --output-path /data/team7_captivevswild/team7_exports

# transferring otu, taxonomy & rooted tree files to local device
scp root@10.19.139.186:/data/team7_captivevswild/team7_exports .
