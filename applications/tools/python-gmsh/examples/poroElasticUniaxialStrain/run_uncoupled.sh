#!/bin/bash

curDir=$PWD
echo $curDir
python3 model_uncoupled.py

rm -rf out-uncoupled
mkdir out-uncoupled
MuPhiSim -input input.inp -inputDir . -outputDir out-uncoupled
