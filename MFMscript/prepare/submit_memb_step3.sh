#!/bin/bash -l
#SBATCH -p panda
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G


singularity run /athena/hwlab/scratch/lab_data/software/apbs_panda/mint17.sif apbs_memb_step3.in boundary.txt >run.log 2>err

exit
