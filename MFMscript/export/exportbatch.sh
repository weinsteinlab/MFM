#!/bin/bash

run_dir="
chol35
chol45
"
run_dir="
cholPS
"
scriptpath="/athena/hwlab/scratch/hex4001/script/MFMscript/export/"
Parentfolder=`pwd`
for i in $run_dir
do
    cd $i
    sbatch $scriptpath/density.sh
    cd $Parentfolder
done

exit

