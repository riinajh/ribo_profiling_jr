#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=filter_undesired_rna
#SBATCH --output=/home/jriina/ribo_profiling_jr/logfiles/08_%j.out
#SBATCH --error=/home/jriina/ribo_profiling_jr/logfiles/08_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

cd /home/jriina/ribo_profiling_jr/raw_data/

module load UHTS/Aligner/bowtie/1.2.0

for file in *_clipped._trimmed.fastq.gz;
do echo $file; gunzip -cd $file | bowtie -S -p 4 ./annotations/undesired_rnas --un $file'filtered.fastq'; #gunzip -c is equal to --stdout, which will output to stdout instead of decompressing the file itself
#-S outputs to sam file, -p specifies number of threads, --un inverts to only nonaligning reads
done
