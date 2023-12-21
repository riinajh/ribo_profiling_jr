#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=filter_undesired_rna
#SBATCH --output=/data/users/jriina/ribo_profiling_jr/logfiles/08_%j.out
#SBATCH --error=/data/users/jriina/ribo_profiling_jr/logfiles/08_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

cd /data/users/jriina/ribo_profiling_jr/raw_data
cp *_clipped_trimmed.fastq.gz ./annotations
cd ./annotations
rm *_clipped._trimmed.fastq.gz

module load UHTS/Aligner/bowtie/1.2.0

for file in $( ls -d *_clipped_trimmed.fastq.gz);
do gunzip $file; done
for file in $( ls -d *_clipped_trimmed.fastq);
do bowtie -S -t -p 4 undesired_rnas ${file} --un ${file:0:24}'_filtered.fastq'; done #gunzip -c is equal to gunzip --stdout, which will output to stdout instead of decompressing the file itself. this is also the same as zcat. 
#-S outputs to sam file, -p specifies number of threads, --un inverts to only nonaligning reads

#for file in $(ls -d *_clipped._trimmed.fastq.gz); do
#gunzip -cd $file | \
#bowtie -p 4 -x undesired_rnas -S  undesired_rnas -t  - --un ${file}_no_r_t_sno_sn_RNA.fastq 2> ${file}_no_r_t_sno_sn_RNA_log.txt > /dev/null;done

#for x in $(ls -d *_clipped_trimmed.fastq.gz); \
#do echo ${x}; \
#gunzip -cd ${x} | \
#bowtie \
#-S \
#-t \
#-p 4 \
#undesired_rnas \
#- \
#--un $(basename ${x} .fastq.gz)_filtered.fastq 2> $(basename ${x} .fastq.gz)_filtered_log.txt > /dev/null; done

#for file in *_clpd_tr.fastq.gz; do
   # gunzip "$file"
#done

#for x in *tr.fastq; do
#    echo ${x}
#    bowtie \
#        -S \
#        -t \
#        -p 4 \
#        ../annotation/undesired_RNA/GRCh38_p13_r_t_sno_sn_RNA_ENSEMBL_NCBI_GtRNAdb \
#        ${x} \
#        --un $(basename ${x} .fastq)_no_r_t_sno_sn_RNA.fastq 2> $(basename ${x} .fastq)_no_r_t_sno_sn_RNA_log.txt > /dev/null
#done

