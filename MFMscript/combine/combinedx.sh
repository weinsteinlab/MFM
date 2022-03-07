#!/bin/bash -l
#SBATCH -p edison
#SBATCH -n1
#SBATCH --mem=20G
conda activate my_opendx

proteinfilename='proteinname'
membranename='membranefoldername'
valency='lipidvalency'
area='lipidarea'
scriptpath='/athena/hwlab/scratch/hex4001/script/MFMscript/combine/'

mkdir $proteinfilename
cp $scriptpath/../run/run_combine/* $proteinfilename/
membenergy=`grep "Total electrostatic energy" ../$membranename/run.log | awk '{ printf "%15.9f\n", $5}'`
systenergy=`grep "Total electrostatic energy" ../protein/run_$proteinfilename.log | awk '{ printf "%15.9f\n", $5}'`
sed -i "s/f_memb_totalenergy/$membenergy/" $proteinfilename/adsorption_F.c
sed -i "s/f_syn_totalenergy/$systenergy/" $proteinfilename/adsorption_F.c
gcc -o $proteinfilename/adsorption_F $proteinfilename/adsorption_F.c -lm
gcc -o $proteinfilename/lipid_mixing $proteinfilename/lipid_mixing.c -lm
gcc -o $proteinfilename/dynamics $proteinfilename/dynamics.c -lm
cp $scriptpath/../run/Irrevelant.pqr .
cp $scriptpath/../run/run_[ad][ny]* $proteinfilename/
cp ../$membranename/boundary.txt $proteinfilename/

cp $scriptpath/combinedx.py combinedx_temp_${proteinfilename}_${membranename}.py
sed -i "s/proteinfilename/$proteinfilename/" combinedx_temp_${proteinfilename}_${membranename}.py
sed -i "s/membranefolder/$membranename/" combinedx_temp_${proteinfilename}_${membranename}.py
python combinedx_temp_${proteinfilename}_${membranename}.py > $proteinfilename/temp 
rm combinedx_temp_${proteinfilename}_${membranename}.py

average_sigma=`grep "average_sigma" $proteinfilename/temp | awk '{ printf "%.15f\n", $3}'`
Debye=`grep "Debye length" ../protein/run_$proteinfilename.log | awk '{ printf "%.5f\n", $3}'`
sed -i "s/average_sigma/$average_sigma/" $proteinfilename/parameters_dynamics.par
sed -i "s/average_sigma/$average_sigma/" $proteinfilename/parameters_mixing.par
sed -i "s/area/$area/" $proteinfilename/parameters_dynamics.par
sed -i "s/area/$area/" $proteinfilename/parameters_mixing.par
sed -i "s/valency/$valency/" $proteinfilename/parameters_dynamics.par
sed -i "s/valency/$valency/" $proteinfilename/parameters_mixing.par
sed -i "s/Debye/$Debye/" $proteinfilename/parameters_mixing.par
rm $proteinfilename/temp

cp $scriptpath/repairmap.sh repairmap_temp_${proteinfilename}_${membranename}.sh
sed -i "s/proteinfilename/$proteinfilename/" repairmap_temp_${proteinfilename}_${membranename}.sh
./repairmap_temp_${proteinfilename}_${membranename}.sh
rm repairmap_temp_${proteinfilename}_${membranename}.sh

exit
