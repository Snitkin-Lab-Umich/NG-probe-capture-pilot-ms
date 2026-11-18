
# Cognac Run
# Start Date: 12/17/24
# Organism: NG
# Activate
library(tidyverse)
library(cognac)
Sys.info()
sessionInfo()

setwd('/nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2024_NG_enrichment_probe_creation/2024-12-17_cognac/results/')

outDir = '/nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2024_NG_enrichment_probe_creation/2024-12-17_cognac/results/'

fnas = readRDS('/nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2024_NG_enrichment_probe_creation/2024-12-17_cognac/data/fna_files.RDS')  
gffs = readRDS('/nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2024_NG_enrichment_probe_creation/2024-12-17_cognac/data/gff_files.RDS')  

fasta_files = paste0('/nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2024_NG_enrichment_probe_creation/2024-12-03_prokka_annotation/results/',fnas)
gff_files = paste0('/nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2024_NG_enrichment_probe_creation/2024-12-03_prokka_annotation/results/',gffs)


parallel::detectCores()

# Run Conac requiring at least 500 genes are included in the alignment. 
algnEnv = cognac(
  fastaFiles = fasta_files ,
  featureFiles = gff_files,
  keepTempFiles = TRUE,
  outDir        = outDir,
  minGeneNum    = 500,
  maxMissGenes  = 0.05,
  threadVal    = 6,
  njTree     = TRUE
)

algnEnv$distMat =CreateAlgnDistMat('/nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2024_NG_enrichment_probe_creation/2024-12-17_cognac/results/concatenated_gene_aa_alignment.fasta', 'raw') 
write.table(algnEnv$distMat, file='/nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2024_NG_enrichment_probe_creation/2024-12-17_cognac/results/distMat.tsv')
 

saveRDS(list(script = "scripts/representative_isolates/02__cognac/02__cognac.R", original_script = "../../2024_NG_enrichment_probe_creation/2024-12-17_cognac/R/cognac.R", generated_at = Sys.time(), inputs = c(fasta_files, gff_files)), file.path(outDir, "provenance.RDS"))
