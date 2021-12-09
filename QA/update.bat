#!/bin/bash
for X in `ls -d TEST*` ; do 
  echo $X 
  cp $X/stdout.* $X/REF/
done


