#!/bin/sh

set -e

# Create new fluid mesh
gmshToFoam fluid.msh -case fluid-openFOAM

# Configure mesh slightly
python3 faces_to_empty.py fluid-openFOAM

cd fluid-openFOAM
rm -rf ?.?
rm -rf ?.??
rm -rf ?.???
touch foam.foam

# create solid input file

cd ../solid-MuPhiSim
rm -rf input output && mkdir input output
cd ../
python3 readFromMesh2D.py
cp input.inp solid-MuPhiSim/input
