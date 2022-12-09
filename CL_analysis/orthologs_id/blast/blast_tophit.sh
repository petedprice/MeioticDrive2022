#!/bin/bash

#$ -l h_rt=8:0:0

#$ -l rmem=32G

#$ -pe smp 8

#$ -P ressexcon
#$ -q ressexcon.q

#$ -wd /fastdata/bop20pp/scRNAseq/wdir

module load apps/python/anaconda2-4.2.0

python --version

python /home/bop20pp/software/MeioticDrive2022/CL_analysis/orthologs_id/blast/21.top-blasthit.py $1 $2
