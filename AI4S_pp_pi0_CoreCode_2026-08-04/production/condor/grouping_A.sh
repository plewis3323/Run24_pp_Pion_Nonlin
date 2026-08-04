#!/bin/bash

# ========== CONFIGURATION ==========
# Define max number of .list files per batch (EDIT THIS VALUE AS NEEDED)
MaxPerBatch=20  # Change this to set the number of .list files per output batch

# Output directory
output_dir="./batched_lists_b"
mkdir -p "$output_dir"

# Find all .list files in the current directory
list_files=($(find . -maxdepth 1 -type f -name "*.list"))
total_files=${#list_files[@]}

# Check if any .list files are found
if [[ "$total_files" -eq 0 ]]; then
    echo "No .list files found in the current directory."
    exit 1
fi

echo "Total .list files found: $total_files"
echo "Processing in batches of $MaxPerBatch files per output list."

batch=1
index=0

# Process in batches of MaxPerBatch
while [[ $index -lt $total_files ]]; do
    batch_file="${output_dir}/output_batch${batch}.list"
    echo "Creating batch file: $batch_file"
    : > "$batch_file"

    # Read up to MaxPerBatch .list files
    for ((i=0; i<MaxPerBatch && index<total_files; i++, index++)); do
        input_file="${list_files[index]}"
        echo "Processing file: $input_file"

        # Read each line in the current .list file
        while IFS= read -r line; do
            # Skip empty lines or comments
            if [[ -z "$line" || "$line" =~ ^# ]]; then
                continue
            fi

            # Extract the last column (assuming it contains the ROOT file)
            root_file=$(echo "$line" | awk '{print $NF}' | xargs)

            # Ensure it is a .root file and format correctly
            if [[ "$root_file" == *.root ]]; then
                echo "$root_file" >> "$batch_file"
            else
                echo "Skipping invalid entry in $input_file: $line"
            fi
        done < "$input_file"
    done

    batch=$((batch + 1))
done

echo "Batch processing complete. Output stored in ${output_dir}/output_batchX.list"

