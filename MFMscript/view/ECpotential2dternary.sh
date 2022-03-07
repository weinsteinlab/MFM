#!/bin/bash -l
#SBATCH -p edison
#SBATCH -n1
#SBATCH --mem=10G

scriptpath="/athena/hwlab/scratch/hex4001/script/MFMscript/view"
conda activate my_opendx
mkdir view

posi='1 500 1000 1600'

for i in $posi
do
gunzip ./$i/cmbChargeDensity.dx.gz
gunzip ./$i/Potential.dx.gz
done

cp $scriptpath/ECpotential2d.py ECpotential2d_temp.py
sed -i "s/posi_list/$posi/" ECpotential2d_temp.py
python ECpotential2d_temp.py
rm ECpotential2d_temp.py

for i in $posi
do
gzip ./$i/cmbChargeDensity.dx
gzip ./$i/Potential.dx
done

exit
