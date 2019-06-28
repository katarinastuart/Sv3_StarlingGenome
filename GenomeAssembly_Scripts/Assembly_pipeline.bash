Assembly-Pipeline--------------------------------------------------[2019.06.13]



https://support.10xgenomics.com/de-novo-assembly/software/overview/latest/performance
https://support.10xgenomics.com/de-novo-assembly/guidance/doc/achieving-success-with-de-novo-assembly
https://support.10xgenomics.com/de-novo-assembly/software/pipelines/latest/using/running


Initial species for sequencing: 

Sturnus vulgaris = HN00105164:

/srv/scratch/z5188231/StarlingGenome/Starling10x/rawdata/HN00105164/HN00105164_10x_RawData_Outs/H2CYFCCX2/fastq_path/H2CYFCCX2/SV01/SV01_S1_L006_R1_001.fastq.gz
/srv/scratch/z5188231/StarlingGenome/Starling10x/rawdata/HN00105164/HN00105164_10x_RawData_Outs/H2CYFCCX2/fastq_path/H2CYFCCX2/SV01/SV01_S1_L006_R2_001.fastq.gz


gunzip -c SV01_S1_L006_R1_001.fastq.gz > /srv/scratch/z5188231/StarlingGenome/Starling10x/data/fastq/SV01_S1_L006_R1_001.fastq
gunzip -c SV01_S1_L006_R2_001.fastq.gz > /srv/scratch/z5188231/StarlingGenome/Starling10x/data/fastq/SV01_S1_L006_R2_001.fastq



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

#run on highmem node

module add supernova/2.1.1

FASTQ=/srv/scratch/z5188231/StarlingGenome/Starling10x/data/fastq/
SAMPLE=SV01

supernova run --id svulgaris-10x --fastqs $FASTQ --sample $SAMPLE --description sturnus_vulgaris_10X --localcores 28 --localmem 900 --maxreads='all'


#I think the fastqs need to be unzipped and sitting in a directory somewhere called fastq, which is what you give to the --fastqs setting. Eg.
supernova run --id brownsnake-supernova2 --fastqs /srv/scratch/babsgenome/BABSGenome-Jun17/raw/2017-09-06.10xHiSeqX/fastq --sample 170814-FR07936203 --description brown_snake_10X-HiSeqX --localcores 28 --localmem 900

supernova run --id rargentea-10x-56X --fastqs $FASTQ --sample $SAMPLE --description rhodamnia_argentea_10x-56X --localcores 28 --localmem 900 --maxreads=450924412

supernova run --id rargentea-10x --fastqs $FASTQ --sample $SAMPLE --description rhodamnia_argentea_10X --localcores 28 --localmem 900 --maxreads='all'