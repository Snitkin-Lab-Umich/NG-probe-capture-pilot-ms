#!/bin/sh
cd /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/data/PathogenWatch_AMR_coverage/database/concatenated_genes
cat /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/data/PathogenWatch_AMR_coverage/database/gene_references/*.fasta > PathogenWatch_reference_genes.fasta

# Make reference dictionary
module load Bioinformatics
module load bwa
bwa index PathogenWatch_reference_genes.fasta