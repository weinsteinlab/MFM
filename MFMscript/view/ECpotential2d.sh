#!/bin/bash -l
#SBATCH -p panda_physbio
#SBATCH -n1
#SBATCH --mem=10G

scriptpath="/athena/hwlab/scratch/hex4001/script/MFMscript/view"
conda activate my_opendx
mkdir view

posi='1 700 1400 2000'

cp $scriptpath/ECpotential2d.py ECpotential2d_temp.py

for i in $posi
do
gunzip ./$i/cmbChargeDensity.dx.gz
gunzip ./$i/Potential.dx.gz
done

sed -i "s/posi_list/$posi/" ECpotential2d_temp.py
python ECpotential2d_temp.py
rm ECpotential2d_temp.py

for i in $posi
do
gzip ./$i/cmbChargeDensity.dx
gzip ./$i/Potential.dx
done

exit
