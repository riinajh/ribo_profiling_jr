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

./Codon_occupancy_cal.sh \
raw_data/annotations/[GRCh38_APPRIS_CDS_18_SingleLine.fa] \
raw_data/annotations/[RPF_KO_Rep1_clpd_tr_no_r_t_sno_sn_RNA_GRCh38_APPRIS_CDS.sam]

mv ./Codon_occupancy.txt ./RPF_KO_Rep1_Codon_occupancy.txt

./Codon_occupancy_cal.sh \
/PATH/reference/GRCh38_APPRIS_CDS_18_SingleLine.fa \
/PATH/RPF_KO_Rep2_clpd_tr_no_r_t_sno_sn_RNA_GRCh38_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_KO_Rep2_Codon_occupancy.txt

./Codon_occupancy_cal.sh \
/PATH/reference/GRCh38_APPRIS_CDS_18_SingleLine.fa \
/PATH/RPF_WT_Rep1_clpd_tr_no_r_t_sno_sn_RNA_GRCh38_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_WT_Rep1_Codon_occupancy.txt

./Codon_occupancy_cal.sh \
/PATH/reference/GRCh38_APPRIS_CDS_18_SingleLine.fa \
/PATH/RPF_WT_Rep2_clpd_tr_no_r_t_sno_sn_RNA_GRCh38_APPRIS_CDS.sam

mv ./Codon_occupancy.txt ./RPF_WT_Rep2_Codon_occupancy.txt
