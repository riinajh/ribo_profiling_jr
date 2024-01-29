#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=genome_mapping
#SBATCH --output=/data/users/jriina/ribo_profiling_jr/logfiles/09_%j.out
#SBATCH --error=/data/users/jriina/ribo_profiling_jr/logfiles/09_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

cd /data/users/jriina/ribo_profiling_jr/raw_data/annotations
module load UHTS/Aligner/bowtie/1.2.0
module load UHTS/Analysis/samtools/1.8
# Now that the reads for the undesired reads have been filtered out, we can align the remaining ones to the genome safely. 
for file in $( ls -d *filtered.fastq); do echo $file; bowtie -S -t -p4 -v1 -m1 --best --strata  annotated_genome -q ${file} | samtools view -h -F4 -b > $(basename ${file} .fastq)_GRCh38.bam; done

for file in $( ls -d *.bam); do echo $file; samtools sort -@ 4 ${file} -o $(basename ${file} .bam)_sorted.bam; done

rm *GRCh38.bam

