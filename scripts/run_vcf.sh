#!/bin/bash

cd /disk/data/vcf/demo
#pwd > $1/log2.txt
singularity exec --bind "/disk/:/disk/" ../development.simg ./vcf2rosetta.py $1$2 >& $1/log.txt
#echo $1$2 >> ${1}/log2.txt
#head -n 5 $1$2 >> ${1}/log2.txtqq
