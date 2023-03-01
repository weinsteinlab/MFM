#!/bin/bash -l
#SBATCH -p panda
#SBATCH -n1
#SBATCH --mem=20G

conda activate my_opendx
python readmap.py

exit
