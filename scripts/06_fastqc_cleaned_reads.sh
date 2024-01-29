#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=qc_cleaned_reads
#SBATCH --output=/home/jriina/ribo_profiling_jr/logfiles/06_%j.out
#SBATCH --error=/home/jriina/ribo_profiling_jr/logfiles/06_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

echo 'Running...'
module add UHTS/Quality_control/fastqc/0.11.7
# Now rerunning fastqc on the clipped+trimmed reads to confirm everything looks better overall.
cd /home/jriina/ribo_profiling_jr/raw_data

fastqc -o /home/jriina/ribo_profiling_jr/raw_data/fastqc -f fastq *_clipped._trimmed.fastq.gz
cd fastqc
rm *.html
for file in {cleaned_RPF_KO_Rep1_fastqc_data.txt,cleaned_RPF_KO_Rep2_fastqc_data.txt,cleaned_RPF_WT_Rep1_fastqc_data.txt,cleaned_RPF_KO_Rep2_fastqc_data.txt}; do touch $file; done
unzip -p  O_Rep1_clipped._trimmed.fastqc.zip O_Rep1_clipped._trimmed_fastqc/fastqc_data.txt > cleaned_RPF_KO_Rep1_fastqc_data.txt
unzip -p  O_Rep2_clipped._trimmed.fastqc.zip O_Rep2_clipped._trimmed_fastqc/fastqc_data.txt > cleaned_RPF_KO_Rep2_fastqc_data.txt
unzip -p  T_Rep1_clipped._trimmed.fastqc.zip T_Rep1_clipped._trimmed_fastqc/fastqc_data.txt > cleaned_RPF_WT_Rep1_fastqc_data.txt
unzip -p  T_Rep2_clipped._trimmed.fastqc.zip T_Rep2_clipped._trimmed_fastqc/fastqc_data.txt > cleaned_RPF_KO_Rep2_fastqc_data.txt
echo 'Finished'
