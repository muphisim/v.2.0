#!/bin/bash

curDir=$PWD
echo $curDir
python3 model.py

Run=1
if [ $Run -eq 1 ]; then
rm -rf $curDir/muphisim  
mkdir $curDir/muphisim
cd $curDir/muphisim && mkdir input output
cd $curDir && cp input.inp $curDir/muphisim/input
cd $curDir/muphisim && MuPhiSim -input input.inp
echo "done simulation in muphisim"
fi

Run=1
if [ $Run -eq 1 ]; then
rm -rf $curDir/muphisim-mpi
mkdir $curDir/muphisim-mpi
cd $curDir/muphisim-mpi && mkdir input output 
cd $curDir && cp input.inp $curDir/muphisim-mpi/input
cd $curDir/muphisim-mpi && mpiexec -np 4 MuPhiSim-mpi -input input.inp
echo "done simulation in muphisim-mpi"
fi
