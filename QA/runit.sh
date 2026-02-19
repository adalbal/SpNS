#!/bin/bash

# Record start time
start_time=$(date +%s)

# Run QA tests
for X in `ls -d TEST_*` ; do
  cd $X
  echo $X

  rm -rf stdout*

  mpirun -np 4 ../../build/a.out params*

  cd ..
done

# Record end time
end_time=$(date +%s)

# Calculate elapsed time
elapsed_time=$((end_time - start_time))

# Print the total execution time
echo -e "\nQA TOTAL EXECUTION TIME: $elapsed_time seconds\n"
