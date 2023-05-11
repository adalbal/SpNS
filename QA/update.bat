#!/bin/bash
for X in `ls -d TEST*` ; do 
  echo $X 
  mkdir -p $X/REF/
  cp $X/stdout.* $X/REF/
done


