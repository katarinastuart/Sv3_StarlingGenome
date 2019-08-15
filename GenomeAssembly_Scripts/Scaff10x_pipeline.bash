#Script for scaffolding short read illumina data with 10x data using Scaff10x
#https://github.com/wtsi-hpag/Scaff10X

    
#!/bin/bash

#PBS -N 2019-08-15.Scaff10x
#PBS -l nodes=1:ppn=32
#PBS -l mem=180gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M katarina.stuart@student.unsw.edu.au
#PBS -m ae

cd /srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.2_Starling10x/scaff/

module add scaff10x

SCAFF10x=/apps/scaff10x/3.1/bin
FASTQ=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.2_Starling10x/data/fastq
READ1=${FASTQ}/SV01_S1_L006_R1_001.fastq
READ2=${FASTQ}/SV01_S1_L006_R2_001.fastq

#barcode extraction
${SCAFF10x}/scaff_BC-reads-1 $READ1 sv10x-BC_1.fastq sv10x-BC_1.name 2>&1 | tee sv10x-BC.log 
${SCAFF10x}/scaff_BC-reads-2 sv10x-BC_1.name $READ2 sv10x-BC_2.fastq 2>&1 | tee -a sv10x-BC.log && touch ~/jobwait/10x-BC.done

#NOTE: Had to fix path to get it to run:
#mkdir scaffdir
#ln -s /share/apps/scaff10x/3.1/bin/* scaffdir
#mkdir scaffdir/scaff-bin
#ln -s /share/apps/scaff10x/3.1/bin/* scaffdir/scaff-bin/

#Run Scaff10x
#Tried a few things for the tiger snake but the basic Minimap mapping without running Break10x seemed to produce the same output as everything else: will just run that for now.

ln /srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.3_GenomeAnnotation/data/GENOME.2018-08-21/Svulgaris_genomic.fna .
SREAD1=sv10x-BC_1.fastq
SREAD2=sv10x-BC_2.fastq
PPN=24
GENOME=Svulgaris_genomic
module add minimap2/2.14
minimap2 -t $PPN -ax sr $PREFIX.fna $SREAD1 $SREAD2 > $PREFIX.10x.sam

CURDIR=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.2_Starling10x/scaff/
${CURDIR}/scaffdir/scaff10x -nodes $PPN -matrix 2000 -reads 8 -score 20 -longread 1 -edge 50000 -link 6 -block 50000 -sam ${CURDIR}/${GENOME}.10x.sam $PREFIX.fna $SREAD1 $SREAD2 $GENOME.scaffolds.fasta 2>&1 | tee sv10x-Scaff10x.log

SCAFF=${GENOME}.scaffolds
python ~/SLiMSuite/tools/seqsuite.py -seqin ${SCAFF}.fasta -log ${SCAFF}.log -summarise -dna i=-1