#!/bin/bash
for X in `ls -d TEST*` ; do
  echo $X
  cd $X
  mkdir -p REF/
  cp stdout.* REF/
  cd ..
done
