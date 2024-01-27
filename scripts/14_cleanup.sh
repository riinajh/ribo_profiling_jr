#!/usr/bin/env bash

#BATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=cleanup
#SBATCH --output=/data/users/jriina/ribo_profiling_jr/logfiles/14_%j.out
#SBATCH --error=/data/users/jriina/ribo_profiling_jr/logfiles/14_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end
for file in `ls -R | grep ".fastq$"`; do gzip $file; done

module load MultiQC/1.11-foss-2021a
cd raw_data/fastqc
multiqc .
cd ../annotations
module load SAMtools/1.13-GCC-10.3.0
for sam in `ls -d *.sam`; do samtools view -h -F 4 -b ${sam} > $(basename ${sam} .fastq)_GRCh38.bam; done

# Sort the BAM file

for x in $(ls -d *.bam); do echo ${x}; samtools sort -@ 4 ${x} -o $(basename ${x} .bam)_sorted.bam; done

# Remove the unsorted BAM file
rm *GRCh38.bam
