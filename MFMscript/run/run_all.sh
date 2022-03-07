#!/bin/bash

run_dir_prefix=""
run_dir="
apo35_1
apo45_1
apoPS_1
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
cholPS_1
cholPS_2
cholPS_3
cholPS_4
cholPS_5
"
start_frame=1
end_frame=1000


scriptpath="/athena/hwlab/scratch/hex4001/script/MFMscript/run/"
Parentfolder=`pwd`
cp $scriptpath/run_dynamics.sh run_dynamics_temp.sh
sed -i "s/START_FRAME/$start_frame/" run_dynamics_temp.sh
sed -i "s/END_FRAME/$end_frame/" run_dynamics_temp.sh
for i in $run_dir
do
    cp run_dynamics_temp.sh ${run_dir_prefix}$i/run_dynamics.sh
    cd ${run_dir_prefix}$i
    sbatch run_dynamics.sh
    cd $Parentfolder
done
rm run_dynamics_temp.sh

exit
