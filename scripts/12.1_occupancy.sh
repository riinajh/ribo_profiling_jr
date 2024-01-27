#! /usr/bin/env bash

#BATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=codon_occupancy
#SBATCH --output=/data/users/jriina/ribo_profiling_jr/logfiles/12_%j.out
#SBATCH --error=/data/users/jriina/ribo_profiling_jr/logfiles/12_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

touch scripts/Codon_occupancy_cal.sh
chmod 774 scripts/Codon_occupancy.sh
wget -O scripts/Codon_occupancy_cal.sh https://raw.githubusercontent.com/LeidelLab/Codon_occupancy_cal/main/Codon_occupancy_cal.sh
