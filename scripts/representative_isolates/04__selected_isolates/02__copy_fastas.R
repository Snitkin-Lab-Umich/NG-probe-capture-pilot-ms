library(tidyverse)

setwd("~/Desktop/ng_mount/")

readRDS("./Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/data/representative_isolates/source/selections/profiles_350.RDS") -> profiles_350
files_350 <- profiles_350$iso_by_qual %>% unname()
files_350 <- paste0(files_350, ".fasta")
file_350_cp_commands <- paste0("cp /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Sequence_data/Public_datasets/NG_fasta_files/",files_350, " /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2024_NG_enrichment_probe_creation/2025-02-24_selected_isolates/results/3000_profiles/fastas/",files_350, " &" )

readRDS("./Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/data/representative_isolates/source/selections/profiles_450.RDS") -> profiles_450
files_450 <- profiles_450$iso_by_qual %>% unname()
files_450 <- paste0(files_450, ".fasta")
file_450_cp_commands <- paste0("cp /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Sequence_data/Public_datasets/NG_fasta_files/",files_450, " /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2024_NG_enrichment_probe_creation/2025-02-24_selected_isolates/results/2500_profiles/fastas/",files_450, " &" )

readRDS("./Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/data/representative_isolates/source/selections/profiles_650.RDS") -> profiles_650
files_650 <- profiles_650$iso_by_qual %>% unname()
files_650 <- paste0(files_650, ".fasta")
file_650_cp_commands <- paste0("cp /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Sequence_data/Public_datasets/NG_fasta_files/",files_650, " /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2024_NG_enrichment_probe_creation/2025-02-24_selected_isolates/results/2000_profiles/fastas/",files_650, " &" )

readRDS("./Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/data/representative_isolates/source/selections/profiles_1050.RDS") -> profiles_1050
files_1050 <- profiles_1050$iso_by_qual %>% unname()
files_1050 <- paste0(files_1050, ".fasta")
files_1050_cp_commands <- paste0("cp /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Sequence_data/Public_datasets/NG_fasta_files/",files_1050, " /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2024_NG_enrichment_probe_creation/2025-02-24_selected_isolates/results/1500_profiles/fastas/",files_1050, " &" )

bin <- "#!/bin/sh"
write_lines(x = c(bin,file_350_cp_commands, file_450_cp_commands, file_650_cp_commands, files_1050_cp_commands),file = "./Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/scripts/representative_isolates/04__selected_isolates/03__cp_bash.sh")
