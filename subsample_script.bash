#!/bin/bash

echo $1 $2 $3 '> echo $1 $2 $3'
#$2 = run number
args=("$@") 

path to subsample script/./subset_fasta.py -f ${args[0]} < 3_mbc_concat.fasta > 3_mbc_concat_${args[0]}%_r${args[1]}.fasta

vsearch --derep_fulllength 3_mbc_concat_${args[0]}%_r${args[1]}.fasta --output 4_mbc_derep_${args[0]}%_r${args[1]}.fasta --sizeout --relabel uniq
rm 3_mbc_subsample_*
