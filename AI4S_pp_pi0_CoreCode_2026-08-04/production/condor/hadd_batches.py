import os
import subprocess

# Parameters
input_dir = "/path/to/root/files"
output_dir = "/path/to/output"
batch_size_limit_gb = 15
batch_prefix = "batch"
max_batch_size_bytes = batch_size_limit_gb * 1024 ** 3

# Get all ROOT files
root_files = []
for root, _, files in os.walk(input_dir):
    for file in files:
        if file.endswith(".root"):
            full_path = os.path.join(root, file)
            root_files.append(full_path)

# Sort files for deterministic batching
root_files.sort()

# Make output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Batching logic
current_batch = []
current_size = 0
batch_index = 0

for file in root_files:
    file_size = os.path.getsize(file)
    if current_size + file_size > max_batch_size_bytes and current_batch:
        # Run hadd on current batch
        output_file = os.path.join(output_dir, f"{batch_prefix}_{batch_index}.root")
        print(f"Merging {len(current_batch)} files into {output_file}...")
        subprocess.run(["hadd", "-f", output_file] + current_batch)
        batch_index += 1
        current_batch = []
        current_size = 0

    current_batch.append(file)
    current_size += file_size

# Handle final batch
if current_batch:
    output_file = os.path.join(output_dir, f"{batch_prefix}_{batch_index}.root")
    print(f"Merging {len(current_batch)} files into {output_file}...")
    subprocess.run(["hadd", "-f", output_file] + current_batch)

