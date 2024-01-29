# In R studio (for example)
setwd('/home/riinajh/HS2023/rna_seq')
# Installation of the package (to be done only once!)
install.packages("devtools")

library("devtools")

install_github(repo = "ohlerlab/RiboseQC")

###### Analysis part ######

# Load the package
library(RiboseQC)

# Prepare genome file (to be done only once!!!)
prepare_annotation_files(annotation_directory = "/home/riinajh/HS2023/rna_seq/qc_analysis/",
                         twobit_file = "/home/riinajh/HS2023/rna_seq/qc_analysis/GRCh38.dna.primary.assembly.2bit",
                         gtf_file = "/home/riinajh/HS2023/rna_seq/qc_analysis/Homo_sapiens.GRCh38.108.gtf",
                         scientific_name = "Homo.sapiens",
                         annotation_name = "GRCh38",
                         export_bed_tables_TxDb = F,
                         forge_BSgenome = T,
                         create_TxDb = T)


# Genome mapped sorted-BAM files

genome_bam <- c("/home/riinajh/HS2023/rna_seq/qc_analysis/WT_Rep1_clipped_trimmed._filtered_GRCh38_sorted_sorted.bam",
                "/home/riinajh/HS2023/rna_seq/qc_analysis/WT_Rep2_clipped_trimmed._filtered_GRCh38_sorted_sorted.bam",
                "/home/riinajh/HS2023/rna_seq/qc_analysis/KO_Rep1_clipped_trimmed._filtered_GRCh38_sorted_sorted.bam",
                "/home/riinajh/HS2023/rna_seq/qc_analysis/KO_Rep2_clipped_trimmed._filtered_GRCh38_sorted_sorted.bam")

load_annotation("/home/riinajh/HS2023/rna_seq/qc_analysis/Homo_sapiens.GRCh38.108.gtf_Rannot")

###### QC plots ######

RiboseQC_analysis(annotation_file ="/home/riinajh/HS2023/rna_seq/qc_analysis/Homo_sapiens.GRCh38.108.gtf_Rannot",
                  bam_files = genome_bam,
                  fast_mode = T,
                  report_file = "/home/riinajh/HS2023/rna_seq/qc_analysis/RPF_samples_QC.html",
                  sample_names = c("WT_Rep1", "WT_Rep2",
                                   "KO_Rep1", "KO_Rep2"),
                  dest_names = c("WT_Rep1", "WT_Rep2",
                                 "KO_Rep1", "KO_Rep2"),
                  write_tmp_files = F)
create_html_report(input_files=c('/home/riinajh/HS2023/rna_seq/WT_Rep1_results_RiboseQC',
                     "/home/riinajh/HS2023/rna_seq/WT_Rep2_results_RiboseQC",
                     "/home/riinajh/HS2023/rna_seq/KO_Rep1_results_RiboseQC",
                     "/home/riinajh/HS2023/rna_seq/KO_Rep2_results_RiboseQC"),
                   input_sample_names=c("WT_R1", "WT_R2",
                                        "KO_R1", "KO_R2"),
                   output_file="/home/riinajh/HS2023/rna_seq/qc_analysis/RPF_samples_QC.html")


