#!/bin/bash -l
#SBATCH -p edison
#SBATCH -n1
#SBATCH --mem=20G

scriptpath="/athena/hwlab/scratch/hex4001/script/MFMscript/view"
conda activate my_opendx
mkdir view

posi='init 1 1500 3500'

for i in $posi
do
gunzip ./$i/cmbChargeDensity.dx.gz
gunzip ./$i/Potential.dx.gz
done

cp $scriptpath/ECpotential.py ECpotential_temp.py
sed -i "s/posi_list/$posi/" ECpotential_temp.py
python ECpotential_temp.py
rm ECpotential_temp.py

for i in $posi
do
gzip ./$i/cmbChargeDensity.dx
gzip ./$i/Potential.dx
done

exit
