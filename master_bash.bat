#!/bin/bash

echo $1 $2 $3 '> echo $1 $2 $3'

args=("$@") 
 
mkdir /home/matthew/test_directory/test_disruns/${args[0]}

cd /home/matthew/test_directory/test_disruns/${args[0]}

echo "Step 1: Denoising" 
vsearch --cluster_unoise /home/matthew/test_directory/testfiles_formaster_bash/4_mbc_derep.fasta --minsize ${args[1]} --unoise_alpha ${args[2]} --centroids 5_mbc_denoise.fasta

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
vsearch --usearch_global home/matthew/files_for_master_bashscript/3_mbc_concat.fasta -db *_otus.fasta -id 0.97 -otutabout ${args[0]}_coleop_reads_map.tsv

echo "Step 8: R Script analysis"
-figure the outputs, which all go into the r script

where does taxonomy come into this? coleoptera asv filtering?

needs all three variables in the r script
-metadata
-reads
-taxonomy

what does the r script output? put into n/n_results?

-also needs adipart script and check_expected_richness

cd../ 
