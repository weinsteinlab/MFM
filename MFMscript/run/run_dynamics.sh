#!/bin/bash -l
#SBATCH -p panda_physbio
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G

#Run
#cd $SGE_O_WORKDIR


for a in `seq START_FRAME END_FRAME`; do

	singularity run /athena/hwlab/scratch/lab_data/software/apbs_panda/mint17.sif apbs.in boundary.txt >run.log 2>err

        ./dynamics parameters_dynamics.par Potential.dx cmbChargeDensity.dx cmbIonAccess.dx ChargeNew.dx

        test=`expr $a % 20`

        if [ $test -eq 0 ] || [ $a -eq 1 ]; then

                mkdir $a

                mv cmbChargeDensity.dx Potential.dx run.log $a

                gzip $a/cmbChargeDensity.dx
                gzip $a/Potential.dx

        fi

        mv ChargeNew.dx cmbChargeDensity.dx

done


