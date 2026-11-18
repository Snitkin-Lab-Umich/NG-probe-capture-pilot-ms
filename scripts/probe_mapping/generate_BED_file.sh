cd /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/data/probe_mapping

bwa index GCF_000006845.1_ASM684v1_genomic.fna 
bwa mem GCF_000006845.1_ASM684v1_genomic.fna 240808-903-Stefan_Green-D10891NGNRH-Ngonorrhoeae.fas > ref-139Kbaits.sam
samtools view -S -b ref-139Kbaits.sam > ref-139Kbaits.bam
bedtools bamtobed < ref-139Kbaits.bam > ref-139Kbaits.bed
sort -k1,1 -k2,2n ref-139Kbaits.bed > ref-139Kbaits.sorted.bed
bedtools merge -i ref-139Kbaits.sorted.bed > ref-139Kbaits.sorted.merge.bed
