#!/bin/bash

echo $1 $2'> echo $1 $2'

args=("$@") 
 
mkdir /home/matthew/diss_runs/${args[0]}

cd /home/matthew/diss_runs/${args[0]}

echo "Step 1: subsample" 
/av/vls/biotools/./subset_fasta.py -f ${args[1]} < /home/matthew/files_mbc_bash/3_mbc_concat.fasta > ${args[0]}_subsample.fasta

echo "Step 2: Dereplication"
vsearch --derep_fulllength ${args[0]}_subsample.fasta --output 4_mbc_derep.fasta --sizeout --sizein --relabel uniq
rm ${args[0]}_subsample.fasta

echo "Step 3: Denoising" 
vsearch --cluster_unoise 4_mbc_derep.fasta --minsize 8 --unoise_alpha 2 --centroids 5_mbc_denoise.fasta
rm 4_mbc_derep.fasta

echo "Step 4: Length Filtering"
vsearch --fastx_filter 5_mbc_denoise.fasta --fastq_minlen 418 --fastq_maxlen 418 -fastaout 6_mbc_indelfil.fasta
rm 5_mbc_denoise.fasta

echo "Step 5: Chimera Filtering"
vsearch --uchime3_denovo 6_mbc_indelfil.fasta --nonchimeras 8_mbc_final.fasta
rm 6_mbc_indelfil.fasta

echo "Step 6: Coleoptera ASV filtering by Thomas's Script"
/av/vls/biotools/./filter_fasta_by_fasta.py k /home/matthew/files_mbc_bash/coleop_asv_filtered.fasta < 8_mbc_final.fasta >9_mbc_beetle_coleop.fasta
rm 8_mbc_final.fasta

echo "Step 7: OTU Delimitation"
vsearch --cluster_size 9_mbc_beetle_coleop.fasta --sizein --relabel otu --id 0.97 --centroids ${args[0]}_otus.fasta
rm 9_mbc_beetle_coleop.fasta

echo "Step 8: Reads Mapping"
vsearch --usearch_global /home/matthew/files_mbc_bash/3_mbc_concat.fasta -db ${args[0]}_otus.fasta -id 0.97 -otutabout ${args[0]}_coleop_reads_map.tsv
mv ${args[0]}_coleop_reads_map.tsv ${args[0]}_reads.tsv
