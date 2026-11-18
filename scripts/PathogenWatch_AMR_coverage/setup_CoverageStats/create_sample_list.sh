path="/scratch/kgontjes_root/kgontjes0/kgontjes/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/data/PathogenWatch_AMR_coverage/CoverageStats/fastq"

sample_id="sample_id"

sample_names=$(ls -1 $path | grep _R1 | cut -d. -f1 | sed 's/_R1_trim_paired//')

echo -e\n $sample_id $sample_names | tr ' ' '\n' > /scratch/kgontjes_root/kgontjes0/kgontjes/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/data/PathogenWatch_AMR_coverage/CoverageStats/config/sample_NG_pilot.tsv

