#!/bin/bash

echo $1 $2

args=("$@") 

mkdir /home/matthew/diss_runs/${args[0]}

cd /home/matthew/diss_runs/${args[0]}

path to subsample script/./subsampling_script subsampleparamater ${args[1]} inputfile 3_mbc_concat outputfile? or...
mv 3_mbc_concat 3_mbc_subsample_${args[0]} 

vsearch --derep_fulllength 3_mbc_subsample_${args[1]}  --output 4_mbc_derep.fasta --sizeout --relabel uniq
rm 3_mbc_subsample_*

#Questions/comments:
#this needs repeating atleast 10 times for each subsampling percentage, therefore I will write a subsampling executing script that 
#calls this script 10 times, similar to the other executing bash script, for each subsample indicating which input file (the output 
#of this one) it will run on, creating a new directory for each call . 
 