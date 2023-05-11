#!/bin/bash

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


