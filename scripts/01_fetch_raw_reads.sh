#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4000M 
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=fetch_raw_reads
#SBATCH --output=/home/jriina/ribo_profiling_jr/logfiles/01_%j.out
#SBATCH --error=/home/jriina/ribo_profiling_jr/logfiles/01_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

if ! [ -e ~/ribo_profiling_jr/raw_data ]; then mkdir ~/ribo_profiling_jr/raw_data; fi
cd /home/jriina/ribo_profiling_jr/raw_data
echo 'now in: '${PWD}
wget -i filelinks.txt #downloading all raw data specified in filelinks.txt

touch samplenames.txt
for name in {RPF_KO_Rep1.fastq.gz,RPF_KO_Rep2.fastq.gz,RPF_WT_Rep1.fastq.gz,RPF_WT_Rep2.fastq.gz}; do echo $name >> samplenames.txt; done
counter=1;for file in ./download.*; do mv $file ./`sed ${counter}'q;d' samplenames.txt`; let counter=counter+1 ; done
#renaming the downloaded files to be more descriptive. a little tricky, this assuming the order of downloaded files stays consistent to the key

