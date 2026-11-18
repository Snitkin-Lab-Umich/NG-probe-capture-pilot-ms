# Use clumpify from bbtools to remove duplicates from fastq files
rule clumpify:
    input:
        r1 = "results/{prefix}/kneaddata/{sample}/{sample}_R1_kneaddata_paired_1.fastq.gz",
        r2 = "results/{prefix}/kneaddata/{sample}/{sample}_R2_kneaddata_paired_2.fastq.gz"
    output:
        r1 = "results/{prefix}/clumpify/{sample}/{sample}_clumpify_R1.fastq.gz",
        r2 = "results/{prefix}/clumpify/{sample}/{sample}_clumpify_R2.fastq.gz"
    benchmark:
        "benchmarks/{prefix}/clumpify/{sample}.benchmark.txt"
    singularity:
         "docker://staphb/bbtools:39.26" 
    threads: 4
    resources:
        mem_mb=8000
    params:
        java_mem="7g"
    shell:
        r"""
        clumpify.sh \
          -Xmx{params.java_mem} \
          t={threads} \
          in1={input.r1:q} \
          in2={input.r2:q} \
          out1={output.r1:q} \
          out2={output.r2:q}
        """