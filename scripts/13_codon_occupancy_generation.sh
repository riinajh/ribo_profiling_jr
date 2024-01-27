#! /usr/bin/env bash

#BATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=codon_occupancy_generation
#SBATCH --output=/data/users/jriina/ribo_profiling_jr/logfiles/13_%j.out
#SBATCH --error=/data/users/jriina/ribo_profiling_jr/logfiles/13_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

cd scripts
echo $PWD
#mv ../raw_data/annotations/isoforms.fa ../raw_data/annotations/isoforms.fa.gz
gunzip ../raw_data/annotations/isoforms.fa.gz

awk '/^>/ { if(NR>1) print "";  printf("%s\n",$0); next; } { printf("%s",$0);}  END {printf("\n");}' < ../raw_data/annotations/isoforms.fa > ../raw_data/annotations/isoforms_single_line.fa

. Codon_occupancy_cal.sh ../raw_data/annotations/isoforms_single_line.fa ../raw_data/annotations/KO_Rep1_clipped_trimmed._filtered_GRCh38_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_KO_Rep1_Codon_occupancy.txt

. Codon_occupancy_cal.sh ../raw_data/annotations/isoforms_single_line.fa ../raw_data/annotations/KO_Rep2_clipped_trimmed._filtered_GRCh38_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_KO_Rep2_Codon_occupancy.txt

. Codon_occupancy_cal.sh ../raw_data/annotations/isoforms_single_line.fa ../raw_data/annotations/WT_Rep1_clipped_trimmed._filtered_GRCh38_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_WT_Rep1_Codon_occupancy.txt

. Codon_occupancy_cal.sh ../raw_data/annotations/isoforms_single_line.fa ../raw_data/annotations/WT_Rep2_clipped_trimmed._filtered_GRCh38_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_WT_Rep2_Codon_occupancy.txt
