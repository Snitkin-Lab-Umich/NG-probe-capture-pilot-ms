# Get to location
cd /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/data/PathogenWatch_AMR_coverage/database

# Download
wget https://raw.githubusercontent.com/pathogenwatch-oss/amr-libraries/refs/heads/main/485.toml

# Parse the toml
conda activate toml_parse     
python3 /nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/scripts/PathogenWatch_AMR_coverage/parse_toml.py