
IsoSeq3-Pipeline--------------------------------------------------[2019.06.13]

https://github.com/PacificBiosciences/IsoSeq3/blob/master/README_v3.1.md

GENBASE=Svulgaris_genomic


#Step 0 - Input
#For each SMRT cell, the movieX.subreads.bam, movieX.subreads.bam.pbi, and movieX.subreadset.xml are needed for processing.



###Step 1 - Circular Consensus Sequence calling
#Each sequencing run is processed by ccs to generate one representative circular consensus sequence (CCS) for each ZMW. Only ZMWs with at least one full pass (at least one subread with SMRT adapter on both ends) are used for the subsequent analysis. Polishing is not necessary in this step and is by default deactivated through.



#!/bin/bash

#PBS -N 2019-06-25.Isoseq1.sh
#PBS -l nodes=1:ppn=24
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

# submitted to katana 
qsub 2019-06-25.Isoseq1

31440.kman.restech.unsw.edu.au



[z5188231@katana2 scripts]$ qsub 2019-06-25.Isoseq1.sh
32685.kman.restech.unsw.edu.au



#For long movies and short inserts, it is advised to limit the number of subreads used per ZMW; this can decrease run-time (only available in ccs version ≥ 3.1.0):
ccs movieX.subreads.bam movieX.ccs.bam --noPolish --minPasses 1 --maxPoaCoverage 10

###Step 2 - Primer removal and demultiplexing
#Removal of primers and identification of barcodes is performed using lima, which offers a specialized --isoseq mode. Even in the case that your sample is not barcoded, primer removal is performed by lima. If there are more than two sequences in your primer.fasta file or better said more than one pair of 5' and 3' primers, please use lima with --peek-guess to remove spurious false positive signal. More information about how to name input primer(+barcode) sequences in this FAQ.
lima movieX.ccs.bam barcoded_primers.fasta movieX.fl.bam --isoseq --no-pbi --peek-guess

###Step 3 - Refine
#Your data now contains full-length reads, but still needs to be refined by:
# - Trimming of poly(A) tails
# - Rapid concatmer identification and removal
#Input The input file for refine is one demultiplexed CCS file with full-length reads and the primer fasta file:
# - <movie.primer--pair>.fl.bam or <movie.primer--pair>.fl.consensusreadset.xml
# - primers.fasta
#Output The following output files of refine contain full-length non-concatemer reads:
# - <movie>.flnc.bam
# - <movie>.flnc.transcriptset.xml

#Actual command to refine:
isoseq3 refine movieX.NEB_5p--NEB_Clontech_3p.fl.bam primers.fasta movieX.flnc.bam

#If your sample has poly(A) tails, use --require-polya. This filters for FL reads that have a poly(A) tail with at least 20 base pairs and removes identified tail:
isoseq3 refine movieX.NEB_5p--NEB_Clontech_3p.fl.bam movieX.flnc.bam --require-polya

###Step 3b - Merge SMRT Cells
#If you used more than one SMRT cells, use dataset for merging. Merge all of your <movie>.flnc.bam files:
dataset create --type TranscriptSet merged.flnc.xml movie1.flnc.bam movie2.flnc.bam movieN.flnc.bam

#Similarly, merge all of your source <movie>.subreadset.xml files:
 dataset create --type SubreadSet merged.subreadset.xml movie1.subreadset.xml movie2.subreadset.xml movieN.subreadset.xml