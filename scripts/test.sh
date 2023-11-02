#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1M
#SBATCH --time=00:01:00
#SBATCH --partition=pall
#SBATCH --job-name=test_script
#SBATCH --output=/home/jriina/ribo_profiling_jr/logfiles/test_%j
#SBATCH --error=/home/jriina/ribo_profiling_jr/logfiles/test_%j
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

echo 'Starting...'
sleep 30
echo '...Ending'
