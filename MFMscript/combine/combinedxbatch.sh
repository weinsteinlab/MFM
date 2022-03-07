#!/bin/bash
filelist="
chol35
chol45
"
filelist="cholPS"
membranename="membrane_PIP2"
membranename="membrane_PS"
lipidarea="70"
lipidarea="65"
lipidvalency="4.0"
lipidvalency="1.0"
scriptpath='/athena/hwlab/scratch/hex4001/script/MFMscript/combine/'
for filename in $filelist
do
  cp $scriptpath/combinedx.sh $scriptpath/combinedx_${filename}.sh
  sed -i "s/proteinname/$filename/" $scriptpath/combinedx_${filename}.sh
  sed -i "s/membranefoldername/$membranename/" $scriptpath/combinedx_${filename}.sh
  sed -i "s/lipidarea/$lipidarea/" $scriptpath/combinedx_${filename}.sh
  sed -i "s/lipidvalency/$lipidvalency/" $scriptpath/combinedx_${filename}.sh
  sbatch $scriptpath/combinedx_${filename}.sh
  rm $scriptpath/combinedx_${filename}.sh
done
