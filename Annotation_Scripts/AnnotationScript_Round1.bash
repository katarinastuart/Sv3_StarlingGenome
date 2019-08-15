GENOME ANNOTATION - S VULGARIS
#test comment

2018-08-21--------------------------------------------------------------[BUSCO-long] ONGOING - RUNNING

# unpack tar.gz file

tar xzvf aves_odb9.tar.gz


#Setting up BUSCO

#!/bin/bash


#PBS -N 2018-08-21.1.BUSCOLong
#PBS -l nodes=1:ppn=16
#PBS -l vmem=124gb
#PBS -l walltime=200:00:00
#PBS -j oe
#PBS -M katarina.stuart@student.unsw.edu.au
#PBS -m ae

GENOMEDIR=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/GENOME.2018-08-21/
GENBASE=Svulgaris_genomic

GENOME=${GENOMEDIR}${GENBASE}.fna
OUTDIR=${GENBASE}.busco
MODEL=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/BUSCO.2018-08-21/aves_odb9

module purge
module add python/3.5.2 blast+/2.2.31 hmmer/3.1b2 augustus/3.2.2 emboss/6.5.7 busco/3.0.2b R/3.3.2-bioconductor-3.5

PATH=$PATH:/share/apps/augustus/3.2.2/scripts/

export BUSCO_CONFIG_FILE="/srv/scratch/z5188231/KStuart.Starling-Aug18/code/BUSCO.config.ini"

mkdir $OUTDIR
cd $OUTDIR


python3 /share/apps/busco/3.0.2b/scripts/run_BUSCO.py -o $GENBASE -i $GENOME -l $MODEL --cpu 16 -m genome --tmp $TMPDIR --long


# submitted job to katana 22/8/18

[z5188231@katana scripts]$ qsub 2018-08-22.1.BUSCOLong.pbs
6679258.katana.science.unsw.edu.au

# timed out after 200 hrs. resubmitting to Katama with the --restart tag

[z5188231@katana scripts]$ qsub 2018-08-22.1.BUSCOLong.pbs
6683522.katana.science.unsw.edu.au


2018-08-22/23--------------------------------------------------------------[Filter transposases] COMPLETE
# Running blast+ and BUSCO config

module add perl/5.20.1
module add maker/2.31.9
module add blast+/2.2.31
module add snap/2013-11-29
module add repeatmasker/4.0.7
module add exonerate/2.2.0
module add python/3.5.2
module add hmmer/3.1b2
module add augustus/3.2.2
module add emboss/6.5.7
module add busco/3.0.2b
module add R/3.3.2-bioconductor-3.5
module add samtools/1.6

# Transposase-protein-database
# Searches the SwissProt transposases.
# Gets the list of SwissProt proteins with matches.
# Generate a list of SwissProt proteins to keep.
# Generate a SwissProt-filtered fasta file.

TpasesPROT=Tpases020812
UniprotSprot=uniprot_sprot.fasta
makeblastdb -in $TpasesPROT -input_type fasta -dbtype prot -out TpasesPROT
blastp -query $UniprotSprot -db TpasesPROT -evalue 1e-10 -max_hsps 1 -max_target_seqs 1 -num_threads 1 -outfmt 6 -out sprot_tpasesprot.tab
cut -f 1 sprot_tpasesprot.tab > sprot_tpaseprot.txt
grep ">" ../TPASES.2018-08-21/uniprot_sprot.fasta | grep -v -f sprot_tpaseprot.txt | sed 's/^>//g' | sed 's/[ ].*//g' > sprot_notpasesprot.txt
xargs samtools faidx $UniprotSprot < sprot_notpasesprot.txt > uniprot_sprot_notpasesprot.fasta


# Transposase-DNA-database
# Now time to do the same with the Transposase DNA database: 

TpasesDNA=Tpases020812DNA
makeblastdb -in $TpasesDNA -input_type fasta -dbtype prot -out TpasesDNA
blastp -query uniprot_sprot_notpasesprot.fasta -db TpasesDNA -evalue 1e-10 -max_hsps 1 -max_target_seqs 1 -num_threads 4 -outfmt 6 -out sprot_tpasesdna.tab
cut -f 1 sprot_tpasesdna.tab > sprot_tpasedna.txt
grep ">" uniprot_sprot_notpasesprot.fasta | grep -v -f sprot_tpasedna.txt | sed 's/^>//g' | sed 's/[ ].*//g' > sprot_clean.txt
xargs samtools faidx uniprot_sprot_notpasesprot.fasta < sprot_clean.txt > uniprot_sprot_clean.fasta

# complete 24/08/18 (took approximately 12-24 hrs to run total)


2018-08-23/24--------------------------------------------------------------[MITES-minature inverted-repeat transposable elements] COMPLETE

# trying to use mitehunter2011
# failed

module add perl/5.20.1
module add mitehunter/2011
module add samtools/1.6

module add muscle
module add mdust
module add formatdb
module add blastall

INPUT=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/GENOME.2018-08-21/Svulgaris_genomic.fna #genome file
PREFIX=Starling
DIR_MITE=/srv/scratch/z5188231/KStuart.Starling-Aug18/programs/MITEhunter/MITEhunter/bin

perl ${DIR_MITE}/MITE_Hunter_manager.pl -i ${INPUT} -g ${PREFIX} - n ${CPU} -S 12345678




perl MITE_Hunter_Installer.pl -f /share/apps/mitehunter/2011/bin/ -d </path/to/dir/with/MITE_Hunter/> -b </path/to/blastall> -M </path/to/muscle> -m </path/to/mdust>


DIR_MITE=/share/apps/mitehunter/2011/bin/

perl /share/apps/mitehunter/2011/bin/MITE_Hunter_manager.pl -i ${INPUT} -g ${PREFIX} -n ${CPU} -S 12345678

perl ${DIR_MITE}/MITE_Hunter_manager.pl -i ${INPUT} -g ${PREFIX} -n ${CPU} -S 12345678
cat *Step8_*.fa > MITE.lib
cp MITE.lib ../


# trying to install and set-up detectmite
# got cd-hit from shared drive

cd /srv/scratch/z5188231/KStuart.Starling-Aug18/programs/detectMITE
tar xzvf cd-hit-v4.6.1-2012-08-27.tgz
mv cd-hit-v4.6.1-2012-08-27 cd-hit
cd cd-hit
make
cd ..
mkdir result

module add matlab/2017b

# first attempt at making detectmite work. Unsuccessful
# incorrect (web) version of cd-hit initially used

data_file=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/GENOME.2018-08-21/Svulgaris_genomic.fna

tic;do_MITE_detection(data_file,'-genome','starling');runtime = toc;

#example script code to input then into matlab through shell - available on the demo file

% Test example for detectMITE.
clear;clc;
data_file = './data/Theobroma_cacao_v1.0.pseudomolecule.fas';
genome_name = 'Theobroma';

tic;
        do_MITE_detection(data_file,'-genome',genome_name,'-cpu',4)
runtime = toc;

fid = fopen('detectMITE.Runtime.txt','a');
fprintf(fid,'-------%s-------\n',genome_name);
fprintf(fid,'Runtime: %f s\n',runtime);
fclose(fid);
exit;


nohup matlab < Test_Demo.m > output.txt &


#starling code for input into matlab
# started running 9:30 am 24/8/18

% Test example for detectMITE.
clear;clc;
data_file = '/srv/scratch/z5188231/KStuart.Starling-Aug18/data/GENOME.2018-08-21/Svulgaris_genomic.fna';
genome_name = 'Starling';

tic;
        do_MITE_detection(data_file,'-genome',genome_name,'-cpu',4)
runtime = toc;

fid = fopen('detectMITE.Runtime.txt','a');
fprintf(fid,'-------%s-------\n',genome_name);
fprintf(fid,'Runtime: %f s\n',runtime);
fclose(fid);
exit;

nohup matlab < Starlingscript.m > Starlingoutput3.txt &

# didn't work. cd-hit-est cannot be found. reshuffling aorund cd-hit folder settup
# let's do this one more time. 24/8/18 at around 10 pm 

nohup matlab < Starlingscript.m > Starlingoutput4.txt &

# i be sad if this doesnt work

# 25/8/18 finished running. Total time of 17.5 hrs. OUTPUTS below
/srv/scratch/z5188231/KStuart.Starling-Aug18/programs/detectMITE
Starling.mite.fasta
Starling.miteSet.fasta
Starlingoutput4.txt

# had to do this because I copied over one of the outputs. same input as 4
nohup matlab < Starlingscript.m > Starlingoutput5.txt &

2018-08-21--------------------------------------------------------------[Maker-with general repeat library] COMPLETE

# Control file settup (maker_opts)

genome=Svulgaris_genomic.fna #genome sequence (fasta file or fasta embeded in GFF3 file)
est=GFDQ01.1.fsa_nt #set of ESTs or assembled mRNA-seq in fasta format
protein=uniprot_sprot.fasta  #protein sequence file in fasta format (i.e. from mutiple oransisms)
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no


# maker script to be sumitted to katana

#!/bin/bash


#PBS -N 2018-08-24.1.makerscript
#PBS -l nodes=38:ppn=4
#PBS -l pvmem=7gb
#PBS -l walltime=12:00:00
#PBS -j oe
#PBS -M katarina.stuart@student.unsw.edu.au
#PBS -m ae

module purge
module load mpich/3.2
module load openmpi/2.0.2
module load perl/5.20.1
module load boost/1.53
module load recon/1.08
module load repeatscout/1.05 
module load trf/4.09
module load rmblast/2.2.28
module load repeatmasker/4.0.7
module load repeatmodeler/1.0.10
module load snap/2013-11-29
module load exonerate/2.2.0
module load genemark-es/4.33
module load trnascan-se/1.3.1
module load blast+/2.6.0
module load maker/2.31.9-mpi

cd /srv/scratch/z5188231/KStuart.Starling-Aug18/data/GENOME.2018-08-21/Svulgaris_genomic.fna

MYGENOME=Starling
OUT_DIR=/srv/scratch/z5188231/KStuart.Starling-Aug18/annotation/2018-08-24.MAKER/

cd ${OUT_DIR}

mpiexec -n 4 maker -base ${MYGENOME} ${OUT_DIR}/maker_opts.ctl ${OUT_DIR}/maker_bopts.ctl ${OUT_DIR}/maker_exe.ctl

# submitted to katana 24/8/18 around 5 pm - subset test run

[z5188231@katana 2018-08-24.MAKER]$ qsub 2018-08-24.makerscript.pbs
6680651.katana.science.unsw.edu.au

#sumbitted to katana 29/8/18 - full run with generic repeat library

[z5188231@katana 2018-08-24.MAKER]$ qsub 2018-08-24.makerscript.pbs
6681865.katana.science.unsw.edu.au

2018-08-30--------------------------------------------------------------[Maker-with aves repeat library] ONGOING - RUNNING

# Change repeat masking model organims to aves... does it work?

model_org=all #select a model organism for RepBase masking in RepeatMasker

# changes to script

cd /srv/scratch/z5188231/KStuart.Starling-Aug18/data/GENOME.2018-08-21/Svulgaris_genomic.fna

MYGENOME=Starling
OUT_DIR=/srv/scratch/z5188231/KStuart.Starling-Aug18/annotation/2018-08-30.MAKER/

cd ${OUT_DIR}

mpiexec -n 4 maker -base ${MYGENOME} ${OUT_DIR}/maker_opts.ctl ${OUT_DIR}/maker_bopts.ctl ${OUT_DIR}/maker_exe.ctl

[z5188231@katana 2018-08-30.MAKER]$ qsub 2018-08-30.makerscript.pbs
6682104.katana.science.unsw.edu.au

# looks like job is on hold as more resources needed for the job. will wait for it to timeout then readdress.
# Resubmitted with new specs

#!/bin/bash


#PBS -N 2018-08-30.1.makerscript
#PBS -l nodes=1:ppn=16
#PBS -l vmem=120gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -M katarina.stuart@student.unsw.edu.au
#PBS -m ae


module purge
module load mpich/3.2
module load openmpi/2.0.2
module load perl/5.20.1
module load boost/1.53
module load recon/1.08
module load repeatscout/1.05
module load trf/4.09
module load rmblast/2.2.28
module load repeatmasker/4.0.7
module load repeatmodeler/1.0.10
module load snap/2013-11-29
module load exonerate/2.2.0
module load genemark-es/4.33
module load trnascan-se/1.3.1
module load blast+/2.6.0
module load maker/2.31.9-mpi

cd /srv/scratch/z5188231/KStuart.Starling-Aug18/data/GENOME.2018-08-21/

MYGENOME=Starling
OUT_DIR=/srv/scratch/z5188231/KStuart.Starling-Aug18/annotation/2018-08-24.MAKER/

cd ${OUT_DIR}

maker -base ${MYGENOME} ${OUT_DIR}/maker_opts.ctl ${OUT_DIR}/maker_bopts.ctl ${OUT_DIR}/maker_exe.ctl

# submitted to katana 4/9/18

[z5188231@katana 2018-08-30.MAKER]$ qsub 2018-08-30.makerscript.pbs
6683867.katana.science.unsw.edu.au




2018-08-27/29--------------------------------------------------------------[LTR (Long Terminal Repeat) retrotransposons] COMPLETE

# unpack tar.gz file

tar xzvf CRL_Scripts1.0.tar.gz 
tar xzvf ProtExcluder1.2.tar.gz

# PART 1------------------------------------------------------ 

module add perl/5.20.1
module add repeatmasker/4.0.7
module add genometools/1.5.9
module add muscle/3.8.31
module add blast+/2.6.0

AR_PATH=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/Results.Part1
GENOME=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/GENOME.2018-08-21/Svulgaris_genomic.fna
PREFIX=Starling
DATE=2018.08.28
DIR_CRL=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/CRL_Scripts1.0
CPU=4
EUK_tRNA=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/eukaryotic-tRNAs.fa
library=${AR_PATH}/LTR/repeats_to_mask_LTR99.fasta

cd ${AR_PATH}
mkdir -p LTR
cd LTR

# File set-up for LTR analysis
awk -F" " '{print $1}' ${GENOME} | sed -e 's/[_.]//g' > ${AR_PATH}/${PREFIX}.${DATE}.nospaces.fasta
gt suffixerator -db ${AR_PATH}/${PREFIX}.${DATE}.nospaces.fasta -indexname ${PREFIX} -tis -suf -lcp -des -ssp -dna

#Find candidate elements with LTRs that are 99%+ similar using LTRharvest
gt ltrharvest -index ${PREFIX} -out ${PREFIX}.out99 -outinner ${PREFIX}.outinner99 -gff3 ${PREFIX}.gff99 -minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 -mintsd 5 -maxtsd 5 -motif tgca -similar 99 -vic 10 > ${PREFIX}.result99

#Find elements with PPT (poly purine tract) or PBS (primer binding site) with LTRdigest
gt gff3 -sort ${PREFIX}.gff99 > ${PREFIX}.gff99.sort
gt ltrdigest -trnas ${EUK_tRNA} ${PREFIX}.gff99.sort ${PREFIX} > ${PREFIX}.gff99.dgt
perl ${DIR_CRL}/CRL_Step1.pl --gff ${PREFIX}.gff99.dgt

​#Additional filtering of the candidate elements
perl ${DIR_CRL}/CRL_Step2.pl --step1 CRL_Step1_Passed_Elements.txt --repeatfile ${PREFIX}.out99 --resultfile ${PREFIX}.result99 --sequencefile ${AR_PATH}/${PREFIX}.${DATE}.nospaces.fasta --removed_repeats CRL_Step2_Passed_Elements.fasta
mkdir -p fasta_files
mv Repeat_*.fasta fasta_files/
mv CRL_Step2_Passed_Elements.fasta fasta_files/
cd fasta_files/
perl ${DIR_CRL}/CRL_Step3.pl --directory ./ --step2 CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25
mv CRL_Step3_Passed_Elements.fasta ../
cd ..


# PART 2------------------------------------------------------ 

module add perl/5.20.1
module add repeatmasker/4.0.7
module add genometools/1.5.9
module add muscle/3.8.31
module add blast+/2.6.0

AR_PATH=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/Results.Part1
PREFIX=Starling
DATE=2018.08.28
INPUT=${AR_PATH}/${PREFIX}.${DATE}.nospaces.fasta
DIR_CRL=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/CRL_Scripts1.0
CPU=4
PBS_NUM_PPN=16
TpasesDNA=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/Tpases020812DNA

# Prepare mite file from earlier analysis
cd /srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/Results.Part1
cp /srv/scratch/z5188231/KStuart.Starling-Aug18/data/MITES.2018.08.27/Starling.miteSet.fasta .
sed '/^-/ d' < Starling.miteSet.fasta > Starling.miteSet.cleaned.fasta
mv Starling.miteSet.cleaned.fasta MITE.lib
cd ${AR_PATH}/LTR

#Identify elements with nested insertions 
perl ${DIR_CRL}/ltr_library.pl --resultfile ${PREFIX}.result99 --step3 CRL_Step3_Passed_Elements.fasta --sequencefile ${AR_PATH}/${PREFIX}.${DATE}.nospaces.fasta
cat lLTR_Only.lib ../MITE.lib > repeats_to_mask_LTR99.fasta
library=${AR_PATH}/LTR/repeats_to_mask_LTR99.fasta
RepeatMasker -pa ${PBS_NUM_PPN} -lib ${library} -nolow -dir . ${AR_PATH}/LTR/${PREFIX}.outinner99
perl ${DIR_CRL}/cleanRM.pl ${PREFIX}.outinner99.out ${PREFIX}.outinner99.masked > ${PREFIX}.outinner99.unmasked
perl ${DIR_CRL}/rmshortinner.pl ${PREFIX}.outinner99.unmasked 50 > ${PREFIX}.outinner99.clean
makeblastdb -in ${TpasesDNA} -dbtype prot
blastx -query ${PREFIX}.outinner99.clean -db ${TpasesDNA} -evalue 1e-10 -num_threads ${CPU} -num_descriptions 10 -out ${PREFIX}.outinner99.clean_blastx.out.txt
perl ${DIR_CRL}/outinner_blastx_parse.pl --blastx ${PREFIX}.outinner99.clean_blastx.out.txt --outinner ${PREFIX}.outinner99


#Building examplars
perl ${DIR_CRL}/CRL_Step4.pl --step3 CRL_Step3_Passed_Elements.fasta --resultfile ${PREFIX}.result99 --innerfile passed_outinner_sequence.fasta --sequencefile ${INPUT}
makeblastdb -in lLTRs_Seq_For_BLAST.fasta -dbtype nucl
blastn -query lLTRs_Seq_For_BLAST.fasta -db lLTRs_Seq_For_BLAST.fasta -evalue 1e-10 -num_descriptions 1000 -out lLTRs_Seq_For_BLAST.fasta.out -num_threads ${CPU}
makeblastdb -in Inner_Seq_For_BLAST.fasta -dbtype nucl
blastn -query Inner_Seq_For_BLAST.fasta -db Inner_Seq_For_BLAST.fasta -evalue 1e-10 -num_descriptions 1000 -out Inner_Seq_For_BLAST.fasta.out -num_threads ${CPU}
perl ${DIR_CRL}/CRL_Step5.pl --LTR_blast lLTRs_Seq_For_BLAST.fasta.out --inner_blast Inner_Seq_For_BLAST.fasta.out --step3 CRL_Step3_Passed_Elements.fasta --final LTR99.lib --pcoverage 90 --pidentity 80



# PART 3------------------------------------------------------

module add perl/5.20.1
module add repeatmasker/4.0.7
module add genometools/1.5.9
module add muscle/3.8.31
module add blast+/2.6.0

AR_PATH=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/Results.Part1
PREFIX=Starling
DATE=2018.08.28
DIR_CRL=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/CRL_Scripts1.0
CPU=4
EUK_tRNA=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/eukaryotic-tRNAs.fa
TpasesDNA=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/Tpases020812DNA
INPUT=${AR_PATH}/${PREFIX}.${DATE}.nospaces.fasta
PBS_NUM_PPN=16

#Old LTRs (85%)
rm fasta_files/* CRL_Step*
gt ltrharvest -index ${PREFIX} -out ${PREFIX}.out85 -outinner ${PREFIX}.outinner85 -gff3 ${PREFIX}.gff85 -minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 -mintsd 5 -maxtsd 5 -vic 10 > ${PREFIX}.result85

#Find elements with PPT or PBS (85%)
gt gff3 -sort ${PREFIX}.gff85 > ${PREFIX}.gff85.sort
gt ltrdigest -trnas ${EUK_tRNA} ${PREFIX}.gff85.sort ${PREFIX} > ${PREFIX}.gff85.dgt
perl ${DIR_CRL}/CRL_Step1.pl --gff ${PREFIX}.gff85.dgt

#Additional filtering of the candidate elements (85%)
perl ${DIR_CRL}/CRL_Step2.pl --step1 CRL_Step1_Passed_Elements.txt --repeatfile ${PREFIX}.out85 --resultfile ${PREFIX}.result85 --sequencefile ${INPUT} --removed_repeats CRL_Step2_Passed_Elements.fasta
mkdir -p fasta_files
mv Repeat_*.fasta fasta_files
mv CRL_Step2_Passed_Elements.fasta fasta_files
cd fasta_files
perl ${DIR_CRL}/CRL_Step3.pl --directory ./ --step2 CRL_Step2_Passed_Elements.fasta --pidentity 60 --seq_c 25
mv CRL_Step3_Passed_Elements.fasta ..
cd ..
#Identify elements with nested insertions (85%)
perl ${DIR_CRL}/ltr_library.pl --resultfile ${PREFIX}.result85 --step3 CRL_Step3_Passed_Elements.fasta --sequencefile ${INPUT}
cat lLTR_Only.lib ../MITE.lib > repeats_to_mask_LTR85.fasta
library=${AR_PATH}/LTR/repeats_to_mask_LTR85.fasta
RepeatMasker -pa ${PBS_NUM_PPN} -lib ${library} -nolow -dir . ${AR_PATH}/LTR/${PREFIX}.outinner85
perl ${DIR_CRL}/cleanRM.pl ${PREFIX}.outinner85.out ${PREFIX}.outinner85.masked > ${PREFIX}.outinner85.unmasked
perl ${DIR_CRL}/rmshortinner.pl ${PREFIX}.outinner85.unmasked 50 > ${PREFIX}.outinner85.clean
blastx -query ${PREFIX}.outinner85.clean -db ${TpasesDNA} -evalue 1e-10 -num_threads ${CPU} -num_descriptions 10 -out ${PREFIX}.outinner85.clean_blastx.out.txt
perl ${DIR_CRL}/outinner_blastx_parse.pl --blastx ${PREFIX}.outinner85.clean_blastx.out.txt --outinner ${PREFIX}.outinner85

#Building examplars (85%).
perl ${DIR_CRL}/CRL_Step4.pl --step3 CRL_Step3_Passed_Elements.fasta --resultfile ${PREFIX}.result85 --innerfile passed_outinner_sequence.fasta --sequencefile ${INPUT}
makeblastdb -in lLTRs_Seq_For_BLAST.fasta -dbtype nucl
blastn -query lLTRs_Seq_For_BLAST.fasta -db lLTRs_Seq_For_BLAST.fasta -evalue 1e-10 -num_descriptions 1000 -out lLTRs_Seq_For_BLAST.fasta.out -num_threads ${CPU}
makeblastdb -in Inner_Seq_For_BLAST.fasta -dbtype nucl
blastn -query Inner_Seq_For_BLAST.fasta -db Inner_Seq_For_BLAST.fasta -evalue 1e-10 -num_descriptions 1000 -out Inner_Seq_For_BLAST.fasta.out -num_threads ${CPU}
perl ${DIR_CRL}/CRL_Step5.pl --LTR_blast lLTRs_Seq_For_BLAST.fasta.out --inner_blast Inner_Seq_For_BLAST.fasta.out --step3 CRL_Step3_Passed_Elements.fasta --final LTR85.lib --pcoverage 90 --pidentity 80

#Consolidate LTRs
<command to shortern LTR names e.g. sed 's/scafNOTSCTS10X2SCAF/TS10X2S/g' ${AR_PATH}/LTR/LTR99.lib | sed 's/dbseq//g' > ${AR_PATH}/LTR/LTR99_edit.lib>
<command to shortern LTR names e.g. sed 's/scafNOTSCTS10X2SCAF/TS10X2S/g' ${AR_PATH}/LTR/LTR85.lib | sed 's/dbseq//g' > ${AR_PATH}/LTR/LTR85_edit.lib>

RepeatMasker -lib LTR99.lib -dir . LTR85.lib
perl ${DIR_CRL}/remove_masked_sequence.pl --masked_elements LTR85.lib.masked --outfile FinalLTR85.lib
cat LTR99.lib FinalLTR85.lib > allLTR.lib


# PART 4------------------------------------------------------

Module purge
module add perl/5.20.1
module load recon/1.08
module load repeatscout/1.05
module load trf/4.09
module load rmblast/2.2.28
module load repeatmasker/4.0.7
module load repeatmodeler/1.0.10

AR_PATH=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/Results.Part1
INPUT=${AR_PATH}/${PREFIX}.${DATE}.nospaces.fasta
PREFIX=Starling
DIR_CRL=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/CRL_Scripts1.0
DIR_PE=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/ProtExcluder1.2
SPROT=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/uniprot_sprot_clean.fasta
<edit path to esl-sfetch e.g. nano mspesl-sfetch.pl>
TpasesPROT=/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/Tpases020812
PBS_NUM_PPN=16



cd $DIR_PE

./Installer.pl -m /share/apps/hmmer/3.1b2/easel/miniapps/ -p /srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/ProtExcluder1.2/

#Repetitive elements with RepeatModeler
cd ${AR_PATH}
cat LTR/allLTR.lib MITE.lib > allMITE_LTR.lib
library=${AR_PATH}/allMITE_LTR.lib
cd ${AR_PATH}/LTR

RepeatMasker -pa ${PBS_NUM_PPN} -lib ${library} -dir . ${INPUT} -q &> repeatmasker.log
perl ${DIR_CRL}/rmaskedpart.pl ${INPUT##*/}.masked 50 > um_${INPUT##*/}
cd ..
mv LTR/um_Starling.2018.08.28.nospaces.fasta .
BuildDatabase -name um_${INPUT##*/}db -engine ncbi um_${INPUT##*/}
nohup RepeatModeler -pa ${PBS_NUM_PPN} -database um_${INPUT##*/}db >& um_${PREFIX}.out

RM_DIR=${AR_PATH}/RM_91438.ThuSep61053272018
perl ${DIR_CRL}/repeatmodeler_parse.pl --fastafile ${RM_DIR}/consensi.fa.classified --unknowns repeatmodeler_unknowns.fasta --identities repeatmodeler_identities.fasta
makeblastdb -in ${TpasesPROT} -dbtype prot
blastx -query repeatmodeler_unknowns.fasta -db ${TpasesPROT} -evalue 1e-10 -num_descriptions 10 -out modelerunknown_blast_results.txt -num_threads ${CPU} 
perl ${DIR_CRL}/transposon_blast_parse.pl --blastx modelerunknown_blast_results.txt --modelerunknown repeatmodeler_unknowns.fasta
mv unknown_elements.txt ModelerUnknown.lib
cat identified_elements.txt repeatmodeler_identities.fasta > ModelerID.lib

#Excluding gene fragments
makeblastdb -in ${SPROT} -dbtype prot
sed 's/(//g' LTR/allLTR.lib | sed 's/)//g' > allLTR_rename.lib
for lib in ModelerID.lib allLTR_rename.lib MITE.lib ModelerUnknown.lib 
do
blastx -query ${lib} -db ${SPROT} -evalue 1e-10 -num_descriptions 10 -num_threads ${CPU} -out ${lib}_blast_results.txt
${DIR_PE}/ProtExcluder.pl ${lib}_blast_results.txt ${lib}
echo -e "${lib}\tbefore\t$(grep -c ">" ${lib})\tafter\t$(grep -c ">" ${lib}noProtFinal)"
done

cat MITE.libnoProtFinal allLTR_rename.libnoProtFinal ModelerID.libnoProtFinal > KnownRepeats.lib
cat KnownRepeats.lib ModelerUnknown.libnoProtFinal > allRepeats.lib


2018-09-04--------------------------------------------------------------[LTR through katana scripts] PAUSED - COMPLETED IN ANOTHER SECTION

# Attempting to run exactly the same way as in notes

/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.09.03

[z5188231@katana RETROTPOSONS.2018.09.03]$ qsub adv_repeats_LTR_part1.pbs
6683866.katana.science.unsw.edu.au

[z5188231@katana RETROTPOSONS.2018.09.03]$ qsub adv_repeats_LTR_part2.pbs
6683869.katana.science.unsw.edu.au

[z5188231@katana RETROTPOSONS.2018.09.03]$ qsub adv_repeats_LTR_part2.pbs
6683874.katana.science.unsw.edu.au

[z5188231@katana RETROTPOSONS.2018.09.03]$ qsub adv_repeats_LTR_part2.pbs
6683878.katana.science.unsw.edu.au

[z5188231@katana RETROTPOSONS.2018.09.03]$ qsub adv_repeats_LTR_part2.pbs
6683880.katana.science.unsw.edu.au


2018-09-07--------------------------------------------------------------[Maker-with species specific repeat library] PAUSED

# ROUND 1------------------------------------------------------

# Copy species specific repeat library from following location
/srv/scratch/z5188231/KStuart.Starling-Aug18/data/RETROTPOSONS.2018.08.27/Results.Part1/allRepeats.lib

# Control file settup (maker_opts)

genome=Svulgaris_genomic.fna #genome sequence (fasta file or fasta embeded in GFF3 file)
est=GFDQ01.1.fsa_nt #set of ESTs or assembled mRNA-seq in fasta format
protein=uniprot_sprot.fasta  #protein sequence file in fasta format (i.e. from mutiple oransisms)
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no

genome=</path/to/genome_assembly.fasta>
protein=</path/to/uniprot_sprot_clean.fasta>
model_org=<e.g. vertebrates or fungi>
rmlib=</path/to/allRepeats.lib>
repeat_protein=</path/to/te_proteins.fasta>
protein2genome=1
trna=1
cpus= 1  #number of cpus requested
min_protein=20
always_complete=1
single_exon=1

# maker script to be sumitted to katana

#!/bin/bash


#PBS -N 2018-08-30.1.makerscript
#PBS -l nodes=1:ppn=16
#PBS -l vmem=120gb
#PBS -l walltime=100:00:00
#PBS -j oe
#PBS -M katarina.stuart@student.unsw.edu.au
#PBS -m ae


module purge
module load mpich/3.2
module load openmpi/2.0.2
module load perl/5.20.1
module load boost/1.53
module load recon/1.08
module load repeatscout/1.05
module load trf/4.09
module load rmblast/2.2.28
module load repeatmasker/4.0.7
module load repeatmodeler/1.0.10
module load snap/2013-11-29
module load exonerate/2.2.0
module load genemark-es/4.33
module load trnascan-se/1.3.1
module load blast+/2.6.0
module load maker/2.31.9-mpi


MYGENOME=Starling
OUT_DIR=/srv/scratch/z5188231/KStuart.Starling-Aug18/annotation/2018-09-07.MAKER/

cd ${OUT_DIR}

maker -base ${MYGENOME} ${OUT_DIR}/maker_opts.ctl ${OUT_DIR}/maker_bopts.ctl ${OUT_DIR}/maker_exe.c$



# submitted to katana 7/9/18 at 10:30am - test run to see if repeat library errors out

[z5188231@katana 2018-09-07.MAKER]$ qsub 2018-09-07.makerscript.pbs
6684629.katana.science.unsw.edu.au
 # now with --restart flag
[z5188231@katana 2018-09-07.MAKER]$ qsub 2018-09-07.makerscript.pbs
6685822.katana.science.unsw.edu.au

# Same thing as above run in /srv/scratch/z5188231/KStuart.Starling-Aug18/annotation/2018-09-09.MAKER/ - FINAL VERSION


2018-09-07--------------------------------------------------------------[Snap Training] ONGOING

module purge
module add python/3.5.2 blast+/2.2.31 hmmer/3.1b2 augustus/3.2.2 emboss/6.5.7 busco/3.0.2b R/3.3.2-bioconductor-3.5

MYGENOME=Starling
MYGENOME_DIR=/srv/scratch/z5188231/KStuart.Starling-Aug18/annotation/2018-09-09.MAKER/
PBS_NUM_PPN=16
LINEAGE=Aves

# ROUND 1------------------------------------------------------

# Create backup of results
tar cvf ${MYGENOME}.maker.output_run1.tar ${MYGENOME}.maker.output/

# Keep gff and fasta files generated by round 1
cd ${MYGENOME_DIR}
mkdir -p results_run1
cd results_run1
gff3_merge -d ../${MYGENOME}.maker.output/${MYGENOME}_master_datastore_index.log
fasta_merge -d ../${MYGENOME}.maker.output/${MYGENOME}_master_datastore_index.log

# Training snap: round 1

cd ${MYGENOME_DIR}/maker/
mkdir -p snap1
cd snap1
ln -s ../results_run1/${MYGENOME}.all.gff ${MYGENOME}.all.gff
maker2zff ${MYGENOME}.all.gff

fathom genome.ann genome.dna -categorize 1000
fathom uni.ann uni.dna -export 1000 -plus
forge export.ann export.dna
hmm-assembler.pl ${MYGENOME} . > ${MYGENOME}.snap1.hmm

# Training Augustus: round 1

python3 /share/apps/busco/3.0.2b/scripts/run_BUSCO.py -o $GENBASE -i $GENOME -l $MODEL --cpu 16 -m genome --tmp $TMPDIR --long

python /share/apps/busco/3.0.2b/scripts/run_BUSCO.py -h

 usage='python BUSCO.py -i [SEQUENCE_FILE] -l [LINEAGE] -o [OUTPUT_NAME] -m [MODE] [OTHER OPTIONS]',

cd ${MYGENOME_DIR}/
run_busco --cpu ${PBS_NUM_PPN} --out ${MYGENOME} --long \
--lineage_path ${LINEAGE} --mode genome \
--in ${MYGENOME}_assembly.fasta

mkdir -p config/species/${MYGENOME}
cp run_${MYGENOME}/augustus_output/retraining_parameters/* config/species/${MYGENOME}/
run_busco --cpu ${PBS_NUM_PPN} --out ${MYGENOME} --long \
--lineage_path ${LINEAGE} --mode genome \
--in ${MYGENOME}_assembly.fasta

# ROUND 2------------------------------------------------------
MYGENOME=Starling
tar cvf ${MYGENOME}.maker.output_run2.tar ${MYGENOME}.maker.output/

mkdir -p results_run2
cd results_run2
gff3_merge -d ../${MYGENOME}.maker.output/${MYGENOME}_master_datastore_index.log
fasta_merge -d ../${MYGENOME}.maker.output/${MYGENOME}_master_datastore_index.log

# At this point it is a good idea to check how many proteins they are predicting
grep -c ">" *.fasta

# Training snap: round 2
mkdir -p snap2
cd snap2
cp ../results_run2/${MYGENOME}.all.gff ./
maker2zff ${MYGENOME}.all.gff

fathom genome.ann genome.dna -categorize 1000
fathom uni.ann uni.dna -export 1000 -plus
forge export.ann export.dna
hmm-assembler.pl ${MYGENOME} . > ${MYGENOME}.snap2.hmm

# CHANGE CONTROL FILES AGAIN
cp maker_opts_run2.ctl maker_opts_run3.ctl
snaphmm=${MYGENOME_DIR}/maker/snap1/${MYGENOME}.snap2.hmm
keep_preds=1

# Final results
cd ${MYGENOME_DIR}/maker/
mkdir -p results_run3
cd results_run3
gff3_merge -d ../${MYGENOME}.maker.output/${MYGENOME}_master_datastore_index.log
fasta_merge -d ../${MYGENOME}.maker.output/${MYGENOME}_master_datastore_index.log


2018-09-24--------------------------------------------------------------[Annotation] ONGOING

module purge
module load mpich/3.2
module load openmpi/2.0.2
module load perl/5.20.1
module load boost/1.53
module load recon/1.08
module load repeatscout/1.05 
module load trf/4.09
module load rmblast/2.2.28
module load repeatmasker/4.0.7
module load repeatmodeler/1.0.10
module load snap/2013-11-29
module load exonerate/2.2.0
module load genemark-es/4.33
module load trnascan-se/1.3.1
module load blast+/2.6.0
module load maker/2.31.9-mpi

MYGENOME=Starling

# Copy files and remove fasta filed from gene preds
mkdir -p annotation
cp * annotation/
cd annotation/
rm *snap* *genemark* *augustus*


# create mapping file
maker_map_ids --prefix SVUL_ --justify 8 ${MYGENOME}.all.gff > ${MYGENOME}.map

# Create *.renamed.fasta and *.renamed.gff files
for i in *.fasta
do
cp ${i} ${i%.fasta}.renamed.fasta
done

cp ${MYGENOME}.all.gff ${MYGENOME}.all.renamed.gff
rm *s.fasta ${MYGENOME}.all.gff #remove the original files in the annotation folder

# rename
map_gff_ids ${MYGENOME}.map ${MYGENOME}.all.renamed.gff

for i in *.renamed.fasta
do
map_fasta_ids ${MYGENOME}.map ${i}
done

# BLAST annotations
makeblastdb -in uniprot_sprot.fasta -input_type fasta -dbtype prot -out uniprot_sprot #get uniprot from master dir

# Split your ${MYGENOME}.all.maker.proteins.renamed.fasta files. This is optional but you can speed this up using a computing cluster and processing in parallel.
mkdir -p split_fasta/
cd split_fasta/
perl fasta-splitter.pl --part-size 1500 --measure count ../${MYGENOME}.all.maker.proteins.renamed.fasta #need to get fasta-splitte.pl from net

BASE_PATH=/srv/scratch/z5188231/KStuart.Starling-Aug18/annotation/2018-09-09.NoBusco.MAKER/results_run3/annotation
FASTA_PATH=${BASE_PATH}/split_fasta
DB=/srv/scratch/z5188231/KStuart.Starling-Aug18/annotation/2018-09-09.NoBusco.MAKER/results_run3/annotation/uniprot_sprot

mkdir -p ${FASTA_PATH}/blast  # BLAST output folder

# As array script
blastp -query ${FASTA_PATH}/${MYGENOME}.all.maker.proteins.renamed.part-${PBS_ARRAYID}.fasta -db ${DB} \
-out ${FASTA_PATH}/blast/${MYGENOME}.all.maker.proteins.renamed.part-${PBS_ARRAYID}.blastout.tsv \
-num_threads 2 -outfmt 6 -evalue 0.000001 -seg yes -soft_masking true -lcase_masking -max_hsps 1


# Merge output of each blast
cd annotation/
cat split_fasta/blast/${MYGENOME}.all.maker.proteins.renamed.part-*.tsv > ${MYGENOME}.all.maker.proteins.renamed.blastout.tsv

# Add the BLAST functional annotation to the

SPROT_FASTA=/srv/scratch/z5188231/KStuart.Starling-Aug18/annotation/2018-09-27.NoBusco.MAKER/results_run3/annotation/uniprot_sprot.fasta

maker_functional_gff ${SPROT_FASTA} ${MYGENOME}.all.maker.proteins.renamed.blastout.tsv ${MYGENOME}.all.renamed.gff > ${MYGENOME}.all.renamed.func.gff
maker_functional_fasta ${SPROT_FASTA} ${MYGENOME}.all.maker.proteins.renamed.blastout.tsv ${MYGENOME}.all.maker.proteins.renamed.fasta > ${MYGENOME}.all.maker.proteins.renamed.func.fasta
maker_functional_fasta ${SPROT_FASTA} ${MYGENOME}.all.maker.proteins.renamed.blastout.tsv ${MYGENOME}.all.maker.transcripts.renamed.fasta > ${MYGENOME}.all.maker.transcripts.renamed.func.fasta
# this line for trnascan line added by me, not in original code?
maker_functional_fasta ${SPROT_FASTA} ${MYGENOME}.all.maker.proteins.renamed.blastout.tsv ${MYGENOME}.all.maker.trnascan.transcripts.renamed.fasta > ${MYGENOME}.all.maker.trnascan.transcripts.renamed.func.fasta


mkdir -p iprs

# InterProScan annotations
# InterProScan is used to add additional protein annotations such as protein families or specific domains (e.g. transmembrane regions). This annotation needs to be performed on the renamed protein fasta file, so we reuse the splitted file.

module load java/8u91
module load perl/5.20.1
module load python/3.6.5
module load signalp/4.1
module load tmhmm/2.0c
module load interproscan/5.25-64.0

BASE_PATH= /srv/scratch/z5188231/KStuart.Starling-Aug18/annotation/2018-09-09.NoBusco.MAKER/annotation
FASTA_PATH=${BASE_PATH}/split_fasta


interproscan.sh -i ${FASTA_PATH}/${MYGENOME}.all.maker.proteins.renamed.part-1.fasta \
-b ${FASTA_PATH}/iprs/${MYGENOME}.all.maker.proteins.renamed.part-1.iprsout \
-cpu 8 -dp -t p -pa -goterms -iprlookup \
-T ${FASTA_PATH}/iprs/tmp \
-appl TIGRFAM,SFLD,Phobius,SUPERFAMILY,PANTHER,Gene3D,Hamap,ProSiteProfiles,Coils,SMART,CDD,PRINTS,ProSitePatterns,SignalP_EUK,Pfam,ProDom,MobiDBLite,PIRSF,TMHMM


# combine outputs
cd results/
cat split_fasta/iprs/${MYGENOME}.all.maker.proteins.renamed.part-*.tsv > ${MYGENOME}.all.maker.proteins.renamed.iprsout.tsv

# We add now the protein domains from InterProScan to the gff file
ipr_update_gff ${MYGENOME}.all.renamed.func.gff ${MYGENOME}.all.maker.proteins.renamed.iprsout.tsv > ${MYGENOME}.all.renamed.func.protdom.gff

# We can also create a track with
iprscan2gff3 ${MYGENOME}.all.maker.proteins.renamed.iprsout.tsv ${MYGENOME}.all.renamed.gff > ${MYGENOME}.all.renamed.visible_domains.gff 




2018-09-03--------------------------------------------------------------[To do stuff] ONGOING - FOR THE REST OF MY LIFE
# see if I can make aves maker run??
# fix the buggy bug bug
# continue running busco long
when busco finishes running - neaten that shit up






