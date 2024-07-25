
This python project is used to couple MuPhiSim & OpenFOAM to run an FSI application. Templates in 2D and 3D are given, with a very similar structure. 

- MuPhiSim: solid simulation
- OpenFOAM: fluid simulation
- gmsh: mesh generation for each coupling iteration
This code additionally requires a new openFOAM library ```libstressTensor.so``` to be coupled. This ```.so``` file needs to be copied to
the directory ```$FOAM_USER_LIBBIN```.

Run ```run.sh``` from terminal to run a coupling simulation. Or ```python3 runFSI.py```

## Configuring the current example
For a given problem ```config.json``` contains parameters used to configure an FSI simulation.

{   
    ## PARAMS THAT CONTROL THE COUPLING
    "caseDirectory": "runFolder",  // folder to run the simulations in.
    "runFromStart": true,  // run a new simulation, or start partway through.
    "maxIterations":10,  // maximum number of iteration the coupling software will run for
    "relativeTolerance": 1e-3,   // convergence tolerance (relative to previous step) 
    "absoluteTolerance": 1e-4  // convergence tolerance (absolute)
    "cleanFiles": true, // removes all but the last iteration files to save space

    ## PARAMS CONTROLLING THE INDIVIDUAL FLUID AND SOLID SIMULATIONS
    "meshSize": 0.15,  // characteristic mesh size for the FSI interface.
    "dimension": 3,  // dimension of the problem.
    "solidProcessors": 6,  // number of processor to run solid simulation.
    "fluidProcessors": [3, 2, 1],  // number of processors to run fluid simulation in each cartesian direction. [num in x dir, num in y dir, num in zdir]. Total number of processors used will be the product of this list. 
    "fluidStressInterfaces": [1],  // list of FSI interfaces between fluid and solid. Indexed by the order they appear in the ```solid.geo``` file.
    "deformingInterfaces": [1, 2],  // list of solid interfaces that deform.
    "youngsModulus": 5e4,  // material properties for the solid
    "poissonRatio": 0.45,  // 
    "solidDensity": 1.05,
    "relaxationModulus": 0.25,  // relaxation factor for solid displacement. 
    "solver": "\"IMPLICIT\"",  // solid solver type: implicit, implicit static. 
    "solverNumSteps": 2000,  // for implicit solvers a set number of steps are given.
}


## Adapting to a new simulation

To create a new problem the mesh files need to be changed. 
User input files include ```fluid.geo``` and ```solid.geo``` which set out the geometry of the fluid and solid demains respectively.
The ```fluid.geo``` file must include the lines

```// Start solid```

```// End solid```

which tells the programme where to input the solid nodes. The last lines of the ```.geo``` files defines the physical regions associated with set of boundary elements in the mesh. In the solid file that order these boundaries are defined gives the index of that boundary.

referenceFluid gives the template directory for each fluid simulation. This must be set up to conform with the standard OpenFOAM directory structure. For 2D problems in OpenFOAM, the mesh is created as 3D with depth 1, and "empty" boundary conditions must be set on the back and front faces of the domain. 

```solidConfig.py``` is the final file that may need user adjustment. Lines specific to a given problem that will need changing are indicated.
