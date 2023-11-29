#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=retry_annotations
#SBATCH --output=/home/jriina/ribo_profiling_jr/logfiles/07_%j.out
#SBATCH --error=/home/jriina/ribo_profiling_jr/logfiles/07_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

cd /home/jriina/ribo_profiling_jr/raw_data/annotations
module load Blast/edirect/2020.08.17

touch undesired.fa
touch tRNAs.fa
wget -O undesired.fa 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query><Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" ><Dataset name = "hsapiens_gene_ensembl" interface = "default" ><Filter name = "transcript_biotype" value = "rRNA,snoRNA,snRNA"/><Attribute name = "ensembl_gene_id" /><Attribute name = "gene_exon_intron" /></Dataset></Query>'
echo 'undesired.fa filesize is currently '`ls -l undesired.fa | cut -c29-39`
wget -O tRNAs.fa http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/hg19-tRNAs.fa
cat tRNAs.fa >> undesired.fa
echo 'undesired.fa filesize is currently '`ls -l undesired.fa | cut -c29-39`
esearch -db nuccore -query "biomol_rRNA [PROP] AND Homo sapiens[ORGN]" | efetch -format fasta >> undesired.fa #needed to append rather than overwrite
echo 'undesired.fa filesize is currently '`ls -l undesired.fa | cut -c29-39`

