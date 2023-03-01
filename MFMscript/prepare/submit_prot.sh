#!/bin/bash -l
#SBATCH -p edison
#SBATCH --cpus-per-task=16
#SBATCH --mem=10G

conda activate apbs_1_5
apbs apbs_prot.in >run_prot.log 2>err_prot

machine=`hostname`

cat <<EOF | sendmail -t
To: hex4001@med.cornell.edu
Subject: Edison job
From: hex4001@med.cornell.edu

Your job on ${machine} has completed.
EOF

exit
