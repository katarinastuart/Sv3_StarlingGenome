#!/bin/bash

#PBS -N 2019-06-25.Isoseq1
#PBS -l nodes=1:ppn=6
#PBS -l mem=7gb
#PBS -l walltime=48:00:00
#PBS -j oe
#PBS -M katarina.stuart@student.unsw.edu.au
#PBS -m ae

#smrtlink isoseq3 

module add smrtlink

DATADIR=/srv/scratch/z5188231/StarlingGenome/StarlingIsoseq/data/20190510_Sequel54261_0015/20190510_Sequel54261_0015/r54261_20190510_034552
OUT_DIR=/srv/scratch/z5188231/StarlingGenome/StarlingIsoseq/analysis

SMRTCELL1=${DATADIR}/1_A01
SMRTCELL2=${DATADIR}/2_B01

cd /srv/scratch/z5188231/StarlingGenome/StarlingIsoseq/scripts/

ccs ${SMRTCELL1}/mm54261_190510_035631.subreads.bam ${OUT_DIR}/m54261_190510_035631.ccs.bam --noPolish --minPasses 1
ccs ${SMRTCELL2}/m54261_190511_001755.subreads.bam ${OUT_DIR}/m54261_190511_001755.ccs.bam --noPolish --minPasses 1

