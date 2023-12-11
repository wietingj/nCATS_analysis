#!/bin/bash

# set partition
#SBATCH -p normal

# set run on x MB node only
#SBATCH --mem 120000

# set run x cpus
#SBATCH --cpus-per-task 24

# set max wallclock time
#SBATCH --time=47:00:00

# set name of job
#SBATCH --job-name=proc

# Add miniconda3 to PATH
. /mnt/frieling/tools/miniconda3/etc/profile.d/conda.sh


# Activate env on cluster node
conda activate nanopolishenv >> /dev/null

echo "Input file: " $1
barcode=$1


# Processing basecalled reads from ONT to generate files for nanopolish.sh 



###############################################################################
# Directory listing
readonly fastq_dir=path/${barcode}
readonly fast5_dir=path
readonly refgenome=/path_to_reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna


#merge all fastq files into one fastq file
    
cat ${fastq_dir}/*fastq > ${barcode}_output.fastq  

#indexing files to raw reads in fast5 format

nanopolish index --directory ${fast5_dir} ${barcode}_output.fastq

#aligning to reference genome using minimap2

minimap2 -L --MD -ax map-ont ${refgenome} ${barcode}_output.fastq 1> ${barcode}_aln.sam

#converting sam to bam

samtools view -S -b ${barcode}_aln.sam 1> ${barcode}_aln.bam

#sorting aln.bam

samtools sort ${barcode}_aln.bam -o ${barcode}_aln.sorted.bam

#indexing
samtools index ${barcode}_aln.sorted.bam
