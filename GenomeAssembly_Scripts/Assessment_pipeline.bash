#Scripts for assembly assessment
#https://github.com/wtsi-hpag/Scaff10X


###Seqsuite--------------------------------------------------[2019.08.19]

    
#!/bin/bash

#PBS -N 2019-08-19.seqsuite
#PBS -l nodes=1:ppn=6
#PBS -l mem=6gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M katarina.stuart@student.unsw.edu.au
#PBS -m ae

cd /srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.2_Starling10x/scaff

GENOME=Svulgaris_genomic
SCAFF=${GENOME}.scaffolds

python ~/SLiMSuite/tools/seqsuite.py -seqin ${SCAFF}.fasta -log ${SCAFF}.log -summarise -dna i=-1



###Busco--------------------------------------------------[2019.08.19]

module purge
module add blast+/2.9.0  python/3.6.5 hmmer/3.2.1 augustus/3.3.2 emboss/6.6.0 busco/3.0.2b

cd /srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.2_Starling10x/scaff

GENOME=Svulgaris_genomic
SCAFF=${GENOME}.scaffolds
BUSCOSET=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.3_GenomeAnnotation/data/BUSCO.2018-08-21
PPN=24

#running busco on just one sample:

#export PATH=/apps/augustus/3.3.2/scripts/:$PATH
python3 /apps/busco/3.0.2b/scripts/run_BUSCO.py -i ${SCAFF}.fasta -o ${SCAFF}.busco -m genome -l ${BUSCOSET}/tetrapoda_odb9/ -c $PPN -sp human

#Running busco on a bunch of samples

#./busco-setup.sh #line needed?

for GENOME in $(ls *.fasta) 
do
    BASEFILE=$(basename $GENOME .fasta)
    echo $GENOME $BASEFILE
    python3 /apps/busco/3.0.2b/scripts/run_BUSCO.py -i ${GENOME} -o $BASEFILE -m genome -l ${BUSCOSET}/tetrapoda_odb9/ -c $PPN -sp human -f
done
