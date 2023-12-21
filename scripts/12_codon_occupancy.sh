#! /usr/bin/env bash

#BATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1000M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=codon_occupancy
#SBATCH --output=/data/users/jriina/ribo_profiling_jr/logfiles/12_%j.out
#SBATCH --error=/data/users/jriina/ribo_profiling_jr/logfiles/12_%j.err
#SBATCH --mail-user=jacob.riina@students.unibe.ch
#SBATCH --mail-type=error,end

#touch scripts/Codon_occupancy_cal.sh
#chmod 774 scripts/Codon_occupancy.sh
#wget -O scripts/Codon_occupancy_cal.sh https://github.com/LeidelLab/Codon_occupancy_cal/blob/main/Codon_occupancy_cal.sh
cd raw_data/annotations

module load UHTS/Aligner/bowtie/1.2.0

for x in $(ls -d *filtered.fastq);do echo ${x};
bowtie -t -p 4 -v 1 -m 1 --best --strata --norc transcriptome -q ${x} -S $(basename ${x} .fastq)_GRCh38_APPRIS_CDS.sam 2> $(basename ${x} .fastq)_GRCh38_APPRIS_CDS_log.txt; done

