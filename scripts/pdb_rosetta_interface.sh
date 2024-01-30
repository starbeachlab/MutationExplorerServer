#!/usr/bin/bash

conda activate pyrose
#cd /scratch/mutationexplorer/rasp-batch
#echo "$(pwd)"
#echo 'python /scratch/mutationexplorer-beta/server/scripts/pdb_rosetta_interface.py'  $1 $2 '-i' $3
python /scratch/mutationexplorer/server/scripts/pdb_rosetta_interface.py $1 $2 -i $3
conda deactivate
