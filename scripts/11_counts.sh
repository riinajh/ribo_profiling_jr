#! /usr/bin/env bash

#BATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=counts_table
#SBATCH --output=/data/users/jriina/ribo_profiling_jr/logfiles/11_%j.out
#SBATCH --error=/data/users/jriina/ribo_profiling_jr/logfiles/11_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

cd raw_data/annotations

module load UHTS/Analysis/subread/2.0.1
module load Subread/2.0.3-GCC-10.3.0
featureCounts -T 8 -t CDS -g gene_id -a Homo_sapiens.GRCh38.108.gtf -o CDS_counts_raw.txt *_sorted_sorted.bam

cut -f 1,7-10 CDS_counts_raw.txt > CDS_counts_processed.txt

featureCounts \
-T 8 \
-t exon \
-g gene_biotype \
-a Homo_sapiens.GRCh38.108.gtf \
-o biotype_counts_rawfile.txt *_sorted_sorted.bam

cut -f 1,7-10 biotype_counts_rawfile.txt > biotype_counts_processed.txt
