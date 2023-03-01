#!/bin/bash -l
#SBATCH -p edison
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G

conda activate apbs_1_5
apbs apbs_memb_step1.in >run_memb_step1.log 2>err_memb_step1

#machine=`hostname`
#
#cat <<EOF | sendmail -t
#To: hex4001@med.cornell.edu
#Subject: Edison job
#From: hex4001@med.cornell.edu
#
#Your job on ${machine} has completed.
#EOF

exit
