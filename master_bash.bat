#!/bin/bash

echo $1 $2 $3 '> echo $1 $2 $3'

args=("$@") 
 
mkdir /home/matthew/diss_runs/${args[0]}

cd /home/matthew/diss_runs/${args[0]}

echo "Step 1: Denoising" 
vsearch --cluster_unoise /home/matthew/files_mbc_bash/4_mbc_derep.fasta --minsize ${args[1]} --unoise_alpha ${args[2]} --centroids 5_mbc_denoise.fasta

echo "Step 2: Length Filtering"
vsearch --fastx_filter 5_mbc_denoise.fasta --fastq_minlen 418 --fastq_maxlen 418 -fastaout 6_mbc_indelfil.fasta
rm 5_mbc_denoise.fasta

echo "Step 3: Translation Filtering"
filtertranslate -i 6_mbc_indelfil.fasta -t 5 
rm *transfail.fa
mv 6_mbc_indelfil_transpass.fa 7_mbc_transpass.fasta
rm 6_mbc_indelfil_transpass.fa

echo "Step 4: Chimera Filtering"
vsearch --uchime3_denovo 7_mbc_transpass.fasta --nonchimeras 8_mbc_final.fasta
rm 7_mbc_transpass.fasta

echo "Step 5: Coleoptera ASV filtering by Thomas's Script"
home/matthew/files_for_master_bashscript./Thomas_filterscript coleoptera_master_asv_file 8_mbc_final.fasta -output 9_mbc_beetle_coleop.fasta 
rm 8_mbc_final.fasta

echo "Step 6: OTU Delimitation"
vsearch --cluster_size 9_mbc_beetle_coleop.fasta --sizein --relabel otu --id 0.97 --centroids ${args[0]}_otus.fasta

echo "Step 7: Reads Mapping"
vsearch --usearch_global /home/matthew/files_mbc_bash/3_mbc_concat.fasta -db *_otus.fasta -id 0.97 -otutabout ${args[0]}_coleop_reads_map.tsv
mv *_coleop_reads_map.tsv reads.tsv
 
# I don't hink the below is 100% correct
echo "Step 8: R Script analysis"
/home/matthew/files_mbc_bash/./analysisscript.R sampledata_matt.csv reads

-also needs adipart script and check_expected_richness

cd../ 
