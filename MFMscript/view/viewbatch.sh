#!/bin/bash

run_dir="
apo35_1
apo45_1
chol35_1
chol35_2
chol35_3
chol35_4
chol35_5
chol45_1
chol45_2
chol45_3
chol45_4
chol45_5
"
run_dir2="
apoPS_1
cholPS_1
cholPS_2
cholPS_3
cholPS_4
cholPS_5
"
scriptpath="/athena/hwlab/scratch/hex4001/script/MFMscript/view/"
Parentfolder=`pwd`
for i in $run_dir
do
    cd $i
    sbatch $scriptpath/ECpotential2d.sh
#    sbatch $scriptpath/density.sh
    cd $Parentfolder
done

exit

