#!/bin/sh

set -e

# Clean up tmeplate file
cd fluid-openFOAM
./clean.sh

rm -rf precice*

cd ../solid-MuPhiSim
rm -rf precice*
rm -rf output/



