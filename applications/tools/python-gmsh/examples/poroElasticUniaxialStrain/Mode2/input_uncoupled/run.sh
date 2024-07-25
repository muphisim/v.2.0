#!/bin/bash

curDir=$PWD
echo $curDir

rm -rf out
mkdir out
MuPhiSim -input input.inp -inputDir . -outputDir out
