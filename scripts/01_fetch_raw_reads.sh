#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=fetch_raw_reads
#SBATCH --output=~/ribo_profiling_jr/logfiles/%j.out
#SBATCH --error=~/ribo_profiling_jr/logfiles/%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

cd ~/ribo_profiling_jr/raw_data
echo $PWD
wget -i filelinks.txt
