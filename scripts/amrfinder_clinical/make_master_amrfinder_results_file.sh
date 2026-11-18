#!/bin/sh

# Usage: script.sh <path_to_results_folder>
# Creates a master AMRFinder TSV from all per-isolate result files in the folder.

if [ -z "$1" ]; then
    echo "Usage: $0 <path_to_results_folder>" fil
    exit 1
fi

results_dir="$1"
filename=$(basename "$results_dir")
master_file="${results_dir}/${filename}_master.tsv"

# Prepend Assembly column to header, skip headers on subsequent files,
# prepend stripped filename to each data row
awk '
    FNR == 1 && NR == 1 {
        print "Assembly\t" $0
        next
    }
    FNR == 1 { next }
    {
        sample = FILENAME
        gsub(/.*\//, "", sample)
        gsub(/\.tsv$/, "", sample)
        print sample "\t" $0
    }
' "$results_dir"/*.tsv > "$master_file"

echo "Created: $master_file"