Assembly-Pipeline--------------------------------------------------[2019.06.13]



https://support.10xgenomics.com/de-novo-assembly/software/overview/latest/performance
https://support.10xgenomics.com/de-novo-assembly/guidance/doc/achieving-success-with-de-novo-assembly
https://support.10xgenomics.com/de-novo-assembly/software/pipelines/latest/using/running


Initial species for sequencing: 

Sturnus vulgaris = HN00105164:

/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.2_Starling10x/rawdata/HN00105164/HN00105164_10x_RawData_Outs/H2CYFCCX2/fastq_path/H2CYFCCX2/SV01/SV01_S1_L006_R1_001.fastq.gz
/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.2_Starling10x/rawdata/HN00105164/HN00105164_10x_RawData_Outs/H2CYFCCX2/fastq_path/H2CYFCCX2/SV01/SV01_S1_L006_R2_001.fastq.gz

gunzip -c SV01_S1_L006_R1_001.fastq.gz > /srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.2_Starling10x/data/fastq/SV01_S1_L006_R1_001.fastq
gunzip -c SV01_S1_L006_R2_001.fastq.gz > /srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.2_Starling10x/data/fastq/SV01_S1_L006_R2_001.fastq



#---
#Ask gunzip to output to standard output and redirect to a file in that directory:
#
#gunzip -c file.gz > /THERE/file
#zcat is a shortcut for gunzip -c.
#
#If you want to gunzip multiple files iterate over all files:
#
#for f in *.gz; do
#  STEM=$(basename "${f}" .gz)
#  gunzip -c "${f}" > /THERE/"${STEM}"
#done
#(here basename is used to get the part of the filename without the extension)
#---


Assembly-supernova--------------------------------------------------[2019.06.28]

#I think the fastqs need to be unzipped and sitting in a directory somewhere called fastq, which is what you give to the --fastqs setting. Eg.
#run on highmem node

module add supernova/2.1.1

FASTQ=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.2_Starling10x/data/fastq/
SAMPLE=SV01

supernova run --id svulgaris-10x --fastqs $FASTQ --sample $SAMPLE --description sturnus_vulgaris_10X --localcores 28 --localmem 900 --maxreads='all'

#[error] The initially computed value for genome size is 1.27 Gb, from which we infer raw coverage of 85.5x. 
#Our recommended range is between 38x and 56x. We test here for extremes of coverage outside the range 30x to 85x. 
#If you have sufficient reads available in the original dataset, we recommend using --maxreads=474082919 (~56x) or a value at least --maxreads=321699123 (~38x).  
#Therefore Supernova is aborting. Your options are: 
#1) Reduce the --maxreads value that you provided on the supernova run  command line, and restart from the beginning. 
#2) Alternatively, if you wish to continue on, without changes, you may do so by restarting supernova run with the option --accept-extreme-coverage, however this is not recommended and not supported. It may cause Supernova to crash, and in most cases will not yield a fully satisfactory assembly. 
#Note that there are reasons why the computed genome size value may not agree with the genome size value that you believe is correct: 
#1) The initially computed value is an approximation. The final value computed by Supernova is likely to be more accurate (although in rare cases it can be far off). 
#2) In extreme cases, episome, microbiome DNA and contamination can significantly affect Supernova's genome size estimate.
#
#Saving pipestance info to svulgaris-10x/svulgaris-10x.mri.tgz
#For assistance, upload this file to 10x Genomics by running:
#
#supernova upload <your_email> svulgaris-10x/svulgaris-10x.mri.tgz

supernova run --id svulgaris-10x-56X --fastqs $FASTQ --sample $SAMPLE --description sturnus_vulgaris_10X-56X --localcores 28 --localmem 900 --maxreads=474082919

module load python/2.7.15

cd /srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.2_Starling10x/assembly/svulgaris-10x-56X/outs
mkdir fasta
supernova mkoutput --asmdir=assembly --outprefix=fasta/svulgaris-10x-470M --style=pseudohap2
cd fasta
gunzip *.gz
python ~/SLiMSuite/tools/seqsuite.py summarise batchrun="*.fasta" basefile=svulgaris-10x-470M dna newlog

#OUTPUT
#File: svulgaris-10x-470M.1.fasta 
#SeqNum: 24417
#TotLength: 1038040562
#MinLength: 1000
#MaxLength: 8802431
#MeanLength: 42513.0
#MedLength: 2458
#N50Length: 695651
#L50Count: 364
#GapLength: 7660150
#GapPC: 0.74
#GCPC: 41.65
#File: svulgaris-10x-470M.2.fasta
#SeqNum: 24417
#TotLength: 1037957698
#MinLength: 1000
#MaxLength: 8799792
#MeanLength: 42509.6
#MedLength: 2458
#N50Length: 695651
#L50Count: 364
#GapLength: 7660120
#GapPC: 0.74
#GCPC: 41.65

#» Predicting a smaller genome size than before » excessive reads » recalculate maxreads and re-run.

#474.11M * 56 / 59.61 = 321.0M » re-run with maxreads=445.39M

FASTQ=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.2_Starling10x/data/fastq/
SAMPLE=SV01

supernova run --id svulgaris-10x-445M --fastqs $FASTQ --sample $SAMPLE --description sturnus_vulgaris_10X-445M --localcores 28 --localmem 900 --maxreads=445390000


cd /srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.2_Starling10x/assembly/svulgaris-10x-445M/outs
mkdir fasta
supernova mkoutput --asmdir=assembly --outprefix=fasta/svulgaris-10x-445M --style=pseudohap2
cd fasta
gunzip *.gz
python ~/SLiMSuite/tools/seqsuite.py summarise batchrun="*.fasta" basefile=svulgaris-10x-445M dna newlog

 Total number of sequences: 26,767
#SUM    00:28:29        Total length of sequences: 1,035,833,441
#SUM    00:28:29        Min. length of sequences: 1,000
#SUM    00:28:29        Max. length of sequences: 7,244,039
#SUM    00:28:29        Mean length of sequences: 38,698.15
#SUM    00:28:29        Median length of sequences: 2,680
#SUM    00:28:29        N50 length of sequences: 504,942
#SUM    00:28:29        L50 count of sequences: 517
#SUM    00:28:29        GC content: 41.63%
#SUM    00:28:29        Gap (N) length: 7,627,620 (0.74%)

 # ~~~~~~ Sequence Summary for svulgaris-10x-445M.2 ~~~~~~~ #
#SUM    00:36:29        Total number of sequences: 447
#SUM    00:36:29        Total length of sequences: 272,603,618
#SUM    00:36:29        Min. length of sequences: 1,574
#SUM    00:36:29        Max. length of sequences: 4,350,710
#SUM    00:36:29        Mean length of sequences: 609,851.49
#SUM    00:36:29        Median length of sequences: 409,361
#SUM    00:36:29        N50 length of sequences: 1,012,626
#SUM    00:36:29        L50 count of sequences: 87
#SUM    00:36:29        GC content: 39.66%
#SUM    00:36:29        Gap (N) length: 1,574,420 (0.58%)
#SAVE   00:36:29        Table "summarise" saved to "svulgaris-10x-445M.summarise.tdt": 2 entries.
#WARN   00:36:29        1 error messages! Check log for details.






####SHARING DIRECTORIES TO STARLING GROUP
chown -R :starling Sv1_StarlingGBS 
chown -R :starling Sv0_GeneralPopgen
chown -R :starling Sv2_StarlingLitReview
chown -R :starling Sv3.2_Starling10x
chown -R :starling Sv3.1_StarlingIsoseq
chown -R :starling Sv3.3_GenomeAnnotation #taking a very long time


chown -R :starling KStuart.Starling-Aug18
ls KStuart.Starling-Aug18 -lrth

chmod 770 KStuart.Starling-Aug18/

ls -lrth | grep Starling

