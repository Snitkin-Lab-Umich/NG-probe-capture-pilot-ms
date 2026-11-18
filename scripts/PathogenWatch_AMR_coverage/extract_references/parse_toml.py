import tomllib
from pathlib import Path

base_dir = Path("/nfs/turbo/umms-esnitkin/Project_NG_surveillance/Analysis/2025_NG_pipeline_testing/NG_probe_capture_ms/data/PathogenWatch_AMR_coverage/database")
output_dir = base_dir / "gene_references"
output_dir.mkdir(exist_ok=True)

with open(base_dir / "485.toml", "rb") as f:
    data = tomllib.load(f)

for gene in data["genes"]:
    with open(output_dir / f'{gene["name"]}.fasta', "w") as out:
        out.write(f'>{gene["name"]}\n{gene["sequence"]}\n')

print(f"Extracted {len(data['genes'])} sequences")