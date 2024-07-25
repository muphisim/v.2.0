#!/bin/bash

curDir=$PWD
echo $curDir


MuPhiSim -input input.inp -FSIConfigFile FSI.config -inputDir input -outputDir output
