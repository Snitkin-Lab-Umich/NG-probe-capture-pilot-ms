directory=/nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/data/NG_probe/kneaddata
output_directory=/scratch/kgontjes_root/kgontjes0/kgontjes/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/data/PathogenWatch_AMR_coverage/CoverageStats/fastq

mkdir -p $output_directory

for sample_dir in $directory/*/; do
    sample=$(basename $sample_dir)
    
    for read in R1 R2; do
        mate=${read#R}
        src="${sample_dir}${sample}_${read}_kneaddata_paired_${mate}.fastq.gz"
        dst="${output_directory}/${sample}_${read}_trim_paired.fastq.gz"
        
        if [ -f "$src" ]; then
            cp "$src" "$dst"
            echo "Copied: $sample $read"
        else
            echo "WARNING: Not found: $src"
        fi
    done
done

echo "Done. Files in output:"
ls $output_directory | wc -l