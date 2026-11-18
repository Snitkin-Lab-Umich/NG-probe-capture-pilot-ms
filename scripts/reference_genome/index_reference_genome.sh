cd /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/data/reference_genome

module load Bioinformatics

module load bwa

bwa index GCF_000006845.1_ASM684v1_genomic.fna

module load samtools

samtools faidx GCF_000006845.1_ASM684v1_genomic.fna

module load gatk

gatk CreateSequenceDictionary -R GCF_000006845.1_ASM684v1_genomic.fna