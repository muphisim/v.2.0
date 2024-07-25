#!/bin/bash

curDir=$PWD
echo $curDir
python3 model.py

rm -rf out
mkdir out
MuPhiSim -input input.inp -inputDir . -outputDir out
