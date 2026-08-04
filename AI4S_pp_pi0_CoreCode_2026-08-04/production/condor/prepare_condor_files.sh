#!/bin/bash

# Create destination directory
mkdir -p Condor_process
> Condor_process/good_runlist.txt

# List of target run numbers
run_numbers=(
54280 54530 54533 54534 54543 54544 54546 54547 54587 54588
54592 54593 54594 54597 54598 54600 54602 54603 54604 54676
54682 54747 54749 54805 54806 54849 54862 54863 54865 54866
54872 54873 54897 54899 54900 54911 54912 54913 54915 54918
54920 54921 54935 54936 54937 54938 54944 54945 54948 54951
54952 54965 54966 54968
)

# Loop through the run numbers
for run in "${run_numbers[@]}"; do
  listfile="dst_calo_run2auau-000${run}.list"

  if [[ -f "$listfile" ]]; then
    cp "$listfile" Condor_process/
    echo "$run" >> Condor_process/good_runlist.txt
  else
    echo "  Missing: $listfile"
  fi
done

echo " Done. Copied available files and created Condor_process/good_runlist.txt"
