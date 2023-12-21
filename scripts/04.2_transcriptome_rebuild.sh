#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000M
#SBATCH --time=03:00:00
#SBATCH --partition=pall
#SBATCH --job-name=transcriptome_rebuild
#SBATCH --output=/data/users/jriina/ribo_profiling_jr/logfiles/04.2_%j.out
#SBATCH --error=/data/users/jriina/ribo_profiling_jr/logfiles/04.2_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

cd /data/users/jriina/ribo_profiling_jr/raw_data/annotations

module load UHTS/Aligner/bowtie/1.2.0

touch transcriptome.fa
wget -O transcriptome.fa 'http://www.ensembl.org/biomart/martservice?query=<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "FASTA" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
			
	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
		<Filter name = "downstream_flank" value = "18"/>
		<Filter name = "upstream_flank" value = "18"/>
		<Filter name = "transcript_biotype" value = "protein_coding"/>
		<Attribute name = "ensembl_transcript_id" />
		<Attribute name = "coding" />
		<Attribute name = "external_gene_name" />
	</Dataset>
</Query>'

bowtie-build transcriptome.fa transcriptome
echo 'transcriptome succeeded'
