#!/bin/bash

for X in `ls -d TEST_*` ; do 
  cd $X
  echo "==============================================================================="
  echo $X 
  echo "==============================================================================="
  cat stdout.0 | grep Final > 1 ; cat REF/stdout.0 | grep Final > 2; diff 1 2; rm -rf 1 2;
  cd ..
done
