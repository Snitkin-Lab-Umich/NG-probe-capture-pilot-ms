
rule bowtie2:
    input:
        r1 = "results/{prefix}/clumpify/{sample}/{sample}_clumpify_R1.fastq.gz",
        r2 = "results/{prefix}/clumpify/{sample}/{sample}_clumpify_R2.fastq.gz",
    output:
        # Alignment to 3k reference set
        bam = temp("results/{prefix}/bowtie2/{sample}/{sample}_sorted.bam"),
        coord_bam = temp("results/{prefix}/bowtie2/{sample}/{sample}_coord.bam"),
        depth = temp("results/{prefix}/bowtie2/{sample}/{sample}.depth"),
        fq1 = "results/{prefix}/bowtie2/{sample}/{sample}_mapped_R1.fastq.gz",
        fq2 = "results/{prefix}/bowtie2/{sample}/{sample}_mapped_R2.fastq.gz",
        full_mapped_reads_info = "results/{prefix}/bowtie2/{sample}/{sample}_bowtie2_full_alignment_info.tsv",

        # Alignment to single-isolate reference
        coverage_bam = temp("results/{prefix}/bowtie2/{sample}/{sample}_coverage_sorted.bam"),
        coverage_coord_bam = temp("results/{prefix}/bowtie2/{sample}/{sample}_coverage_coord.bam"),
        coverage_depth = temp("results/{prefix}/bowtie2/{sample}/{sample}_coverage.depth"),
        coverage_mapped_reads_info = "results/{prefix}/bowtie2/{sample}/{sample}_coverage_alignment_info.tsv",
        cov_row = temp("results/{prefix}/bowtie2/{sample}/{sample}_coverage_row.tsv")
    params:
        index = config["bowtie2_index"],
        index_coverage = config["bowtie2_index_for_coverage"],
        sample = "{sample}"
    benchmark:
        "benchmarks/{prefix}/bowtie2/{sample}.benchmark.txt"
    conda:
        "envs/bowtie2_samtools.yaml"
    shell:
        """
        # Map original reads to the 3k reference set
        bowtie2 -x {params.index} \
            -1 {input.r1} \
            -2 {input.r2} \
            -p 10 \
        | samtools view -Sb -F 4 - \
        | samtools sort -n -@ 6 -o {output.bam} -

        samtools sort \
            -o {output.coord_bam} \
            {output.bam}

        samtools index {output.coord_bam}

        # Depth from the 3k reference alignment
        samtools depth \
            -aa \
            -q 20 \
            -Q 20 \
            {output.coord_bam} \
        > {output.depth}

        # Mapping statistics for the 3k reference alignment
        samtools flagstat \
            {output.bam} \
            -O tsv \
        > {output.full_mapped_reads_info}

        # Convert mapped reads to paired end FASTQs
        samtools fastq \
            -1 >(gzip > {output.fq1}) \
            -2 >(gzip > {output.fq2}) \
            -0 /dev/null \
            -s /dev/null \
            -n \
            {output.bam}

        # Map extracted FASTQs to the single isolate reference
        bowtie2 -x {params.index_coverage} \
            -1 {output.fq1} \
            -2 {output.fq2} \
            -p 10 \
        | samtools view -Sb -F 4 - \
        | samtools sort -n -@ 6 -o {output.coverage_bam} -

        samtools sort \
            -o {output.coverage_coord_bam} \
            {output.coverage_bam}

        samtools index {output.coverage_coord_bam}

        # Calculate depth from the single isolate alignment
        samtools depth \
            -aa \
            -q 20 \
            -Q 20 \
            {output.coverage_coord_bam} \
        > {output.coverage_depth}

        # Mapping statistics for the single-isolate alignment
        samtools flagstat \
            {output.coverage_bam} \
            -O tsv \
        > {output.coverage_mapped_reads_info}

        # Calculate coverage using the single-isolate reference length
        GENOME_LEN=$(awk '{{s+=$2}} END{{print s}}' \
            {params.index_coverage}.fasta.fai)

        echo -e "Sample\t>=1x\t>=5x\t>=10x\t>=25x\t>=50x" \
            > {output.cov_row}

        ROW="{params.sample}"

        for T in 1 5 10 25 50; do
            COV_BASES=$(awk -v T=$T \
                '$3>=T{{c++}} END{{print c+0}}' \
                {output.coverage_depth})

            COVERAGE=$(awk \
                -v c=$COV_BASES \
                -v L=$GENOME_LEN \
                'BEGIN{{printf "%.2f", (c/L*100)}}')

            ROW="$ROW\t$COVERAGE"
        done

        echo -e "$ROW" >> {output.cov_row}
        """


rule merge_cov_summary:
    input:
        bowtie2_rows = expand("results/{prefix}/bowtie2/{sample}/{sample}_coverage_row.tsv", prefix=PREFIX, sample=SAMPLES),
    output:
        bowtie2_summary = "results/{prefix}/{prefix}_Report/bowtie2_coverage_summary.tsv",
    shell:
        """
        # Write static header
        echo -e "Sample\t>=1x\t>=5x\t>=10x\t>=25x\t>=50x" > {output.bowtie2_summary}
        
        # Append all per-sample rows (skip headers in each row file)
        tail -n +2 -q {input.bowtie2_rows} >> {output.bowtie2_summary}
        """
