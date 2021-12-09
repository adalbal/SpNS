#!/bin/bash
for X in `ls -d TEST_*` ; do 
  cd $X
  echo $X 
  rm -rf OUTPUT/*
  rm -rf slurm* stdout* output*
  sbatch runit.cmd
  cd ..
done


