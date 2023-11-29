#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=clip_trim_raw_reads
#SBATCH --output=/home/jriina/ribo_profiling_jr/logfiles/05_%j.out
#SBATCH --error=/home/jriina/ribo_profiling_jr/logfiles/05_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

cd /home/jriina/ribo_profiling_jr/raw_data

module load UHTS/Quality_control/cutadapt/2.5

for file in *.fastq.gz; do cutadapt -j 4 -q 25 -m 25 -a CTGTAGGCACCATCAAT --trimmed-only -o ${file:4:7}'_clipped.fastq.gz' ${file}; done
for file in *_clipped.fastq.gz; do cutadapt -j 4 -q 25 -m 25 -u -3 -o ${file:1:15}'_trimmed.fastq.gz' ${file}; done



