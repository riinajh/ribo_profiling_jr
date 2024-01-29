#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=qc_raw_reads
#SBATCH --output=/home/jriina/ribo_profiling_jr/logfiles/02_%j.out
#SBATCH --error=/home/jriina/ribo_profiling_jr/logfiles/02_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

echo 'Running...' #performing inital fasqc analysis on the raw reads
module add UHTS/Quality_control/fastqc/0.11.7
mkdir /home/jriina/ribo_profiling_jr/raw_data/fastqc
cd /home/jriina/ribo_profiling_jr/raw_data
fastqc -o /home/jriina/ribo_profiling_jr/raw_data/fastqc -f fastq *.fastq.gz
cd fastqc
rm *.html
for file in (RPF_KO_Rep1_fastqc_data.txt,RPF_KO_Rep2_fastqc_data.txt,RPF_WT_Rep1_fastqc_data.txt,RPF_KO_Rep2_fastqc_data.txt); do touch $file; done
#for file in *; unzip -p $file `echo $file | cut -c1-18`'/fastqc_data.txt' > `echo $file | cut -c1-18`'_data.txt'
unzip -p  RPF_KO_Rep1_fastqc.zip RPF_KO_Rep1_fastqc/fastqc_data.txt > RPF_KO_Rep1_fastqc_data.txt
unzip -p  RPF_KO_Rep2_fastqc.zip RPF_KO_Rep2_fastqc/fastqc_data.txt > RPF_KO_Rep2_fastqc_data.txt
unzip -p  RPF_WT_Rep1_fastqc.zip RPF_WT_Rep1_fastqc/fastqc_data.txt > RPF_WT_Rep1_fastqc_data.txt
unzip -p  RPF_WT_Rep2_fastqc.zip RPF_WT_Rep2_fastqc/fastqc_data.txt > RPF_KO_Rep2_fastqc_data.txt
echo 'Finished'
#specifies output dir, filetype, then files to operate over. option to generate extracted files instead of zips
