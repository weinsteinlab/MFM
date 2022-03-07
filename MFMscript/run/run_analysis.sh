#!/bin/bash -l
#SBATCH -p panda_physbio
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G

for dir in `seq START_FRAME 20 END_FRAME`; do

gunzip $dir/cmbChargeDensity.dx.gz
grep Total $dir/run.log | awk '{ printf "%15.9f\n", $5}' >> el_energy_apbs.dat
./lipid_mixing parameters_mixing.par $dir/cmbChargeDensity.dx >> lipid_mixing.dat
gzip $dir/cmbChargeDensity.dx

done

wc -l el_energy_apbs.dat > temp
awk '{ printf "%d\n", $1 }' temp > parameters_af.dat

awk '{ printf $11 " " $18 "\n"}' lipid_mixing.dat > test.dat
./adsorption_F parameters_af.dat test.dat el_energy_apbs.dat > adsorption_F.dat
rm test.dat
rm temp
rm parameters_af.dat
#gnuplot plot_energy.dat


#./check_pdb combinedcharge.dx parameters_dynamics.dat pip2_fraction2.pdb analysis2.pdb > pip2_number2.dat

exit
