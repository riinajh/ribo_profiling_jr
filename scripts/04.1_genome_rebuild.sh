#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=03:00:00
#SBATCH --partition=pall
#SBATCH --job-name=genome_rebuild
#SBATCH --output=/data/users/jriina/ribo_profiling_jr/logfiles/04.1_%j.out
#SBATCH --error=/data/users/jriina/ribo_profiling_jr/logfiles/04.1_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

cd /data/users/jriina/ribo_profiling_jr/raw_data/annotations

module load UHTS/Aligner/bowtie/1.2.0

gunzip genome_primary_assembly.fa.gz
bowtie-build genome_primary_assembly.fa annotated_genome
echo 'genome succeeded'
