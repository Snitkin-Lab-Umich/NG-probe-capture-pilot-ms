import os
import math

# Constants
total_isolates = 33055
num_jobs = 509
isolates_per_job = total_isolates // num_jobs  # 65 isolates per job
output_dir = "/nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/scripts/representative_isolates/01__prokka_annotation/sbatch_jobs"

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Model script as a template
model_script = """#!/bin/sh
# Job name
#SBATCH --job-name={job_name}
# User info
#SBATCH --mail-user=aryansin@umich.edu
#SBATCH --mail-type=BEGIN,END
#SBATCH --export=ALL
#SBATCH --partition=standard
#SBATCH --account=esnitkin1
#SBATCH --nodes=1  --ntasks=1 --cpus-per-task=10 --mem=50g --time=12:00:00

IFS=$'\\n'       # make newlines the only separator
set -f          # disable globbing

# Isolates of interest
isolates=`ls /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Sequence_data/Public_datasets/NG_fasta_files | cut -d. -f1 | head -n {end} | tail -n {count}`

# Initiate a list of commands
command_list_array=()

# Run for loop, where it generates isolate command for each of the forward end files.

for isolate in $isolates
do
        # Echo isolate name
        echo $isolate
        # Create current command
        curr_cmd="prokka --outdir /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2024_NG_enrichment_probe_creation/2024-12-03_prokka_annotation/results/$isolate --prefix $isolate  --addgenes --locustag $isolate --kingdom Bacteria --genus Neisseria --species gonorrhoeae  --centre X --compliant /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Sequence_data/Public_datasets/NG_fasta_files/$isolate.fasta "
        echo $curr_cmd
        # Add the prokka command to list
        command_list_array+=("$curr_cmd")
done

# Run convict commands using the unix parallel command.
parallel ::: "${{command_list_array[@]}}"
"""

# Generate sbatch jobs
for job_num in range(1, num_jobs + 1):
    start = (job_num - 1) * isolates_per_job + 1
    end = start + isolates_per_job - 1
    
    # Handle the last job (in case isolates don't divide evenly)
    if job_num == num_jobs:
        end = total_isolates
    
    job_name = f"job_{job_num}"
    script_content = model_script.format(job_name=job_name, end=end, count=isolates_per_job)
    output_file = os.path.join(output_dir, f"{job_name}.sbat")
    
    # Write the script to a file
    with open(output_file, "w") as f:
        f.write(script_content)

print(f"Generated {num_jobs} sbat jobs in {output_dir}")