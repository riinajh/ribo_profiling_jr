#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=download_annotations
#SBATCH --output=/home/jriina/ribo_profiling_jr/logfiles/03_%j.out
#SBATCH --error=/home/jriina/ribo_profiling_jr/logfiles/03_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

mkdir /home/jriina/ribo_profiling_jr/raw_data/annotations
cd /home/jriina/ribo_profiling_jr/raw_data/annotations
module load Blast/edirect/2020.08.17

touch undesired.fa.gz #downloading all fasta sequences for our undeisred rnas(we filter out later)
touch tRNAs.fa
wget -O undesired.fa.gz http://www.ensembl.org/biomart/martresults/54?file=martquery_1103150659_772.txt.gz
gunzip undesired.fa.gz
wget -O tRNAs.fa http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/hg19-tRNAs.fa
cat tRNAs.fa >> undesired.fa
esearch -db nuccore -query "biomol_rRNA[prop] AND "Homo sapiens"[Organism]" | efetch -format fasta > undesired.fa

touch genome_primary_assembly.fa.gz # and now the human genome and annotations
touch Homo_sapiens.GRCh38.108.gtf.gz
wget -O genome_primary_assembly.fa.gz https://ftp.ensembl.org/pub/release-108/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -O Homo_sapiens.GRCh38.108.gtf.gz https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz
gunzip genome_primary_assembly.fa.gz

touch isoforms.fa # finally, transcripts
wget -O isoforms.fa http://www.ensembl.org/biomart/martresults/55?file=martquery_1113165952_922.txt.gz
