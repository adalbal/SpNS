#!/bin/bash

TESTNUMBER=$1

if [ -z "$TESTNUMBER" ]; then
  echo "Usage: $0 <test_number>"
  exit 1
fi

# Check if test exists
if [ -z "$(ls -d TEST_${TESTNUMBER}_* 2>/dev/null)" ]; then
  echo "Error: No test TEST_${TESTNUMBER}_* found"
  exit 1
fi

for X in `ls -d TEST_${TESTNUMBER}_*` ; do
  echo "Updating: $X"
  cd $X
  mkdir -p REF/
  cp stdout.* REF/
  cd ..
done
