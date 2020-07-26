#!/bin/bash

echo $1 $2 $3 '> echo $1 $2 $3'

args=("$@") 
 
mkdir /home/matthew/diss_runs/${args[0]}

cd /home/matthew/diss_runs/${args[0]}
if [ ${args[2]} == X ]; then
"Step 1: No Alpha, therefore minsize only"
vsearch --derep_fulllength /home/matthew/files_mbc_bash/4_mbc_derep.fasta --output 5_mbc_denoise.fasta --sizeout --relabel uniq --minuniquesize ${args[1]}
else 
echo "Step 1: Denoising" 
vsearch --cluster_unoise /home/matthew/files_mbc_bash/4_mbc_derep.fasta --minsize ${args[1]} --unoise_alpha ${args[2]} --centroids 5_mbc_denoise.fasta
fi

echo "Step 2: Length Filtering"
vsearch --fastx_filter 5_mbc_denoise.fasta --fastq_minlen 418 --fastq_maxlen 418 -fastaout 6_mbc_indelfil.fasta
rm 5_mbc_denoise.fasta

echo "Step 3: Chimera Filtering"
vsearch --uchime3_denovo 7_mbc_transpass.fasta --nonchimeras 8_mbc_final.fasta
rm 7_mbc_transpass.fasta

echo "Step 4: Coleoptera ASV filtering by Thomas's Script"
home/matthew/files_for_master_bashscript/./filter_fasta_by_fasta.py k coleoptera_master_asvs.fasta < 8_mbc_final.fasta >9_mbc_beetle_coleop.fasta
rm 8_mbc_final.fasta

echo "Step 5: OTU Delimitation"
vsearch --cluster_size 9_mbc_beetle_coleop.fasta --sizein --relabel otu --id 0.97 --centroids ${args[0]}_otus.fasta
rm 9_mbc_beetle_coleop.fasta

echo "Step 6: Reads Mapping"
vsearch --usearch_global /home/matthew/files_mbc_bash/3_mbc_concat.fasta -db ${args[0]}_otus.fasta -id 0.97 -otutabout ${args[0]}_coleop_reads_map.tsv
mv ${args[0]}_coleop_reads_map.tsv reads.tsv

