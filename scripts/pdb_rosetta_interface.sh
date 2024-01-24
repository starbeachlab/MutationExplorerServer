#!/usr/bin/bash

echo "################ interface.sh before conda #######################"
conda activate pyrose
#cd /scratch/mutationexplorer/rasp-batch
echo "######################### interface.sh ##################################"
python /scratch/mutationexplorer/server/scripts/pdb_rosetta_interface.py  $1 $2 -i $3
conda deactivate
