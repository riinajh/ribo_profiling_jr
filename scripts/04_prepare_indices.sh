#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=03:00:00
#SBATCH --partition=pall
#SBATCH --job-name=prepare_indices
#SBATCH --output=/home/jriina/ribo_profiling_jr/logfiles/04_%j.out
#SBATCH --error=/home/jriina/ribo_profiling_jr/logfiles/04_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

cd /home/jriina/ribo_profiling_jr/raw_data/annotations

module load UHTS/Aligner/bowtie/1.2.0

bowtie-build undesired.fa undesired_rnas
echo 'undesired rnas succeeded'
bowtie-build genome_primary_assembly.fa.gz annotated_genome
echo 'genome succeeded'
bowtie-build isoforms.fa transcript_isoforms
echo 'isoforms succeeded'

