#!/bin/bash

# Record start time
start_time=$(date +%s)

# Run QA tests
for X in `ls -d TEST_*` ; do 
  cd $X
  echo $X 

  rm -rf OUTPUT/*
  rm -rf slurm* stdout* output*

  if which sbatch >/dev/null; then
    sbatch runit.cmd
  else
    bash runit.cmd
  fi

  cd ..
done

# Record end time
end_time=$(date +%s)

# Calculate elapsed time
elapsed_time=$((end_time - start_time))

# Print the total execution time
echo -e "\nQA TOTAL EXECUTION TIME: $elapsed_time seconds\n\n"
