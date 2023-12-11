#!/bin/bash
## This script is for nanopolish on slurm

# set partition
#SBATCH -p normal

# set run on x MB node only
#SBATCH --mem 120000

# set run x cpus
#SBATCH --cpus-per-task 24

# set max wallclock time
#SBATCH --time=47:00:00

# set name of job
#SBATCH --job-name=nanop

# Add miniconda3 to PATH
. /mnt/frieling/tools/miniconda3/etc/profile.d/conda.sh


###Script for calling methylation using nanopolish

# Activate environment where nanopolish is installed
conda activate nanopolishenv >> /dev/null

########Define Input, Gene Names and Positions
echo "input argument: " $1
barcode=$1

####### Input of Basecalling/
input_fastq=path/${barcode}_output.fastq    ####output.fastq
input_bam=path/${barcode}_aln.sorted.bam    ###aln.sorted.bam

refgenome=/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

###########
echo "test sleep"


nanopolish call-methylation \
    -t 32 \
    -r ${input_fastq} \
    -b ${input_bam} \
    -g ${refgenome} \
    1> ${barcode}_methylation_calls.tsv
