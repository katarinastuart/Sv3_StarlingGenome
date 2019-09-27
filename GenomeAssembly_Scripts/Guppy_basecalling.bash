###Guppy basecalling


#https://community.nanoporetech.com/protocols/Guppy-protocol/v/gpb_2003_v1_revm_14dec2018/run-guppy-on-linux
#With config file:
#  guppy_basecaller -i <input path> -s <save path> -c <config file> [options]
#flow cell R9.4.1
#on GPU add: -num_callers 8 -ipc_threads 16


DATA=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.3_StarlingNanopore/data/fast5
OUTPUT=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.3_StarlingNanopore/data/basecall

module load guppy/3.2.1

guppy_basecaller --print_workflows
guppy_basecaller -h

# -i input file
# --qscore_filtering
#
#
#

#guppy_basecaller -i $DATA -s $OUTPUT -c dna_r9.4.1_450bps_hac.cfg --gpu_runners_per_device 4 --num_callers 8 --compress_fastq --qscore_filtering --min_qscore 7 


guppy_basecaller -i $DATA -s $OUTPUT -c dna_r9.4.1_450bps_hac.cfg --cpu_threads_per_caller 12 --num_callers 8 --compress_fastq --qscore_filtering --min_qscore 7 

#Started running at 3:41

# --resume

###Guppy basecalling


#https://community.nanoporetech.com/protocols/Guppy-protocol/v/gpb_2003_v1_revm_14dec2018/run-guppy-on-linux
#With config file:
#  guppy_basecaller -i <input path> -s <save path> -c <config file> [options]
#flow cell R9.4.1
#on GPU add: -num_callers 8 -ipc_threads 16


DATA=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.3_StarlingNanopore/data/fast5.2
OUTPUT=/srv/scratch/z5188231/KStuart.Starling-Aug18/Sv3.3_StarlingNanopore/data/basecall.2

module load guppy/3.2.1

guppy_basecaller --print_workflows
guppy_basecaller -h

# -i input file
# --qscore_filtering
#
#
#

#guppy_basecaller -i $DATA -s $OUTPUT -c dna_r9.4.1_450bps_hac.cfg --gpu_runners_per_device 4 --num_callers 8 --compress_fastq --qscore_filtering --min_qscore 7 


guppy_basecaller -i $DATA -s $OUTPUT -c dna_r9.4.1_450bps_hac.cfg --cpu_threads_per_caller 12 --num_callers 8 --compress_fastq --qscore_filtering --min_qscore 7 

#Started running at 3:41

# --resume