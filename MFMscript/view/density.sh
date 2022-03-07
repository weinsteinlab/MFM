#!/bin/bash -l
#SBATCH -p panda_physbio
#SBATCH -n1
#SBATCH --mem=20G

scriptpath="/athena/hwlab/scratch/hex4001/script/MFMscript/view"
conda activate my_opendx
mkdir view

posi='1 700 1400 2000'

cp $scriptpath/density.py density_temp.py

for i in $posi
do
gunzip ./$i/cmbChargeDensity.dx.gz
done

sed -i "s/posi_list/$posi/" density_temp.py
python density_temp.py
rm density_temp.py

for i in $posi
do
gzip ./$i/cmbChargeDensity.dx
done

exit
