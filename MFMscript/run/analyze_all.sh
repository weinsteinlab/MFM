#!/bin/bash

run_dir_prefix="apo_"
#run_dir="f l r u"
#run_dir="fu fd fr fl fru frd flu fld"
run_dir="fu fru flu" #PIP2
run_dir="fu fl fru flu fld" #PIP
#run_dir="fu fl fru flu fld" #PS
#run_dir_prefix="fast_dist_"
#run_dir="3 4 6 8 10 15"
start_frame=1020
end_frame=1500


scriptpath="/athena/hwlab/scratch/hex4001/script/MFMscript/run/"
Parentfolder=`pwd`
cp $scriptpath/run_analysis.sh run_analysis_temp.sh
sed -i "s/START_FRAME/$start_frame/" run_analysis_temp.sh
sed -i "s/END_FRAME/$end_frame/" run_analysis_temp.sh
for i in $run_dir
do
    cp run_analysis_temp.sh ${run_dir_prefix}$i/run_analysis.sh
    cd ${run_dir_prefix}$i
    sbatch run_analysis.sh
    cd $Parentfolder
done
rm run_analysis_temp.sh

exit
