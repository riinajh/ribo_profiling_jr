#! /usr/bin/env bash

#BATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=2bit
#SBATCH --output=/data/users/jriina/ribo_profiling_jr/logfiles/10_%j.out
#SBATCH --error=/data/users/jriina/ribo_profiling_jr/logfiles/10_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end
cd raw_data/annotations
module load SequenceAnalysis/blat/36

faToTwoBit genome_primary_assembly.fa GRCh38.dna.primary.assembly.2bit
