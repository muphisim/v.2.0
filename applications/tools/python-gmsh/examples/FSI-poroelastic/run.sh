#!/bin/sh

set -e

# Clean up tmeplate file
cd fluid-openFOAM
./clean.sh

cd ../
./config.sh

cd fluid-openFOAM
./run.sh


