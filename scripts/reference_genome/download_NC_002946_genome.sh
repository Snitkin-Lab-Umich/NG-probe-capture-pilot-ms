conda activate ncbi_datasets_v_18_33_1
cd /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/data/reference_genome
datasets download genome accession GCF_000006845.1 --filename NG_reference_genome.zip
unzip NG_reference_genome.zip
mv /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/data/reference_genome/ncbi_dataset/data/GCF_000006845.1/GCF_000006845.1_ASM684v1_genomic.fna GCF_000006845.1_ASM684v1_genomic.fna
rm -rf ncbi_dataset