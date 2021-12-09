#!/bin/bash

for X in `ls -d TEST_*` ; do 
  cd $X
  echo "==============================================================================="
  echo $X 
  echo "==============================================================================="
  cat stdout.0 | grep Time > 1 ; cat REF/stdout.0 | grep Time > 2; diff 1 2; rm -rf Time 1 2;
  cd ..
done
