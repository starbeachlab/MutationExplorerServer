#!/bin/bash

cd /scratch/mutationexplorer/vcf/demo
#pwd > $1/log2.txt
singularity exec --bind "/scratch/:/scratch/" ../development.simg ./vcf2rosetta.py $1$2 >& $1/log.txt
#echo $1$2 >> ${1}/log2.txt
#head -n 5 $1$2 >> ${1}/log2.txtqq
