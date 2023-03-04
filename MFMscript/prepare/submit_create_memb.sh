#!/bin/bash -l
#SBATCH -p edison
#SBATCH -n1
#SBATCH --mem=20G
conda activate my_opendx

python create_memb.py 
