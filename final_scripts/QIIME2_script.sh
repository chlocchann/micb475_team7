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

# removed “/“ from metadata column names and replaced with “_” in Excel, then exported  as zoo_metadata.tsv
# upload zoo metadata.tsv file into working directory
scp zoo_metadata.tsv root@10.19.139.186:/data/team7_captivevswild 


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

# create a detached screen session & activate QIIME2 to run overnight
screen -S zoo_dada2
conda activate qiime2-2023.7

# determining ASVs with DADA2
# based on qiime2 visualization of demux.qzv, truncation length of 233bp was selected
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
  --m-sample-metadata-file /data/team7_captivevswild/zoo_metadata.tsv
  
qiime feature-table tabulate-seqs \
  --i-data rep_seqs.qza \
  --o-visualization rep_seqs.qzv

# transferring visualization files to local device
scp root@10.19.139.186:/data/team7_captivevswild/denoise_stats.qzv .
scp root@10.19.139.186:/data/team7_captivevswild/zoo_table.qzv .
scp root@10.19.139.186:/data/team7_captivevswild/rep_seqs.qzv .

## TRAINING & USING CLASSIFIER FOR TAXONOMIC ANALYSIS ##

# create a detached screen session, extract V4 rRNA seqs from Silva & run overnight
screen -S zoo_silva
conda activate qiime2-2023.7

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
  --i-reads rep_seqs.qza \
  --o-classification taxonomy.qza

# visualizing taxonomy file
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# transferring taxonomy file to local device
scp root@10.19.139.186:/data/team7_captivevswild/taxonomy.qzv .

## FILTERING ##

# filtering out mitochondrial and chloroplast samples
qiime taxa filter-table \
  --i-table zoo_table.qza \
  --i-taxonomy taxonomy.qza \
  --p-exclude mitochondria,chloroplast \
  --o-filtered-table zoo_filtered_table.qza

# filtering out features with 5 or lower counts
qiime feature-table filter-features \
  --i-table zoo_filtered_table.qza \
  --p-min-frequency 5 \
  --o-filtered-table zoo_frequency_filtered_table.qza

qiime feature-table summarize \
  --i-table zoo_frequency_filtered_table.qza \
  --o-visualization zoo_frequency_filtered_table.qzv \
  --m-sample-metadata-file /data/team7_captivevswild/zoo_metadata.tsv

# transferring filtered table visualization files to local device
scp root@10.19.139.186:/data/team7_captivevswild/zoo_frequency_filtered_table.qzv .

## PHYLOGENETIC DIVERSITY ANALYSES ##

# create a detached screen session
screen -S zoo_phylo
conda activate qiime2-2023.7

# generating a tree for phylogenetic diversity analyses
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep_seqs.qza \
  --o-alignment aligned-rep_seqs.qza \
  --o-masked-alignment masked-aligned-rep_seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza 

# alpha-rarefaction
qiime diversity alpha-rarefaction \
  --i-table zoo_filtered_table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 214000 \
  --m-metadata-file /data/team7_captivevswild/zoo_metadata.tsv \
  --o-visualization alpha-rarefaction.qzv

# transferring alpha rarefaction file to local device
scp root@10.19.139.186:/data/team7_captivevswild/alpha-rarefaction.qzv .

# calculating alpha- & beta-diversity metrics
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table zoo_filtered_table.qza \
  --p-sampling-depth 83290 \
  --m-metadata-file /data/team7_captivevswild/zoo_metadata.tsv \
  --output-dir core-metrics-results

## PICRUSt2 ANALYSIS ##

# create a detached screen session
screen -S zoo_picrust
conda activate qiime2-2023.7

# running the picrust2 qiime2 plugin
qiime picrust2 full-pipeline \
  --i-table zoo_frequency_filtered_table.qza \
  --i-seq rep_seqs.qza \
  --output-dir picrust2_output \
  --p-placement-tool sepp \
  --p-hsp-method pic \
  --p-max-nsti 2 \
  --verbose

## EXPORTING OTU, TAXONOMY, ROOTED TREE & PICRUSt2 PATHWAY FILES ##

# creating directory for export
mkdir team7_exports

# OTU table export
qiime tools export \
  --input-path /data/team7_captivevswild/zoo_filtered_table.qza \
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

# picrust2 pathway export
qiime tools export \
  --input-path picrust2_output/pathway_abundance.qza \
  --output-path /data/team7_captivevswild/team7_exports/picrust2_output

biom convert -i feature-table.biom --to-tsv -o pathway_abundance.tsv

# transferring otu, taxonomy & rooted tree files to local device
scp root@10.19.139.186:/data/team7_captivevswild/team7_exports .
