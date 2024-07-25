import os
full_path = os.path.realpath(__file__)
print(full_path + "\n")

import gmsh
import GmshModel
import Configure
import CreateInputFile



class myConfigure (Configure.Configure):
    def __init__(self):
        """
        configure in oxfem
        Attributes:
            FEM     : all FEM locations
            MM      : all MM locations
            solvers : all solvers
            materials: all material behaviors
                
        """
        # define FEM domains
        self.FEM = {
                   #support All or [(dim,phys), ...] or nothing
                   "Support": "All"
                   #"Support": [(3,1)]
                   }
        # define MM domains
        self.MM = {
              #support All or [(dim,phys), ...] or nothing
              #"Support": "All"
              #"Support": [(3,1)]
              }
        # define solvers - a tuple of dictionary
        self.solvers = ({
                #solver type
                "Solver": "IMPLICIT STATIC",
                # number of step
                "Number Of Steps": 1000,
                #scale factor to estimate time sep for solver
                "Relative Tolerance": 1e-8,
                "FSI": ((1,1),),
                #output specification, each component of the list corresponds a line in the oxfemm input file
                "Outputs": ((1000,), ("All" , "max"),),
                # start and end time
                "Time": (0.,100.),
                "DispBC":(
                            # each bc contains (type, location(dim, physical) list of prescried components and prescribed values
                            ("DISPLACEMENT INST",(1,4), ((0,0.),(1,0.))), 
                            ("DISPLACEMENT INST",(1,2), ((1,0.),)), 
                            ("DISPLACEMENT INST",(1,3), ((1,0.),)), 
#                            ("DISPLACEMENT INST",(1,1), ((1,0.),)), 
                         ), 
                "Extraction":(
                               ("NODE","Displacement_x","Unknown","Rough",(1,1),0),
                               #("NODE","CenterDisplacement","Unknown","Rough",(1,16),3),
                            )
                },
                )
        # define all Neumann BCs
        self.NeumannBCs = {
                "PressureBC": (
                                # type, location, time interval, value
#                                 ("PRESSURE RAMP", (1,1), (0.,5), -0.4),
#                                 ("PRESSURE INST", (1,1), (5.,10.), -0.4),
                              ),
                "Flux ExtraDof": (
                                    # type, location, time interval, value
                                    #("HEAT INST", (2,52), (0., 0.5), (800.,800.,0.)),
                                    # type, location, time interval, value
                                    #("VOLUMETRIC HEAT FLUX", (3,1), (0.,1.), 100.),
                                    # type, location, time interval, h and sink temperature
                                    #("CONVECTION", (2,52), (0., 0.5), (800.,0)),
                                    # type, location, time interval, A and sink temperature
                                    #("RADIATION", (2,52), (0., 0.5), (500.,0)),
                                  )
                #"GRAVITY": 1.
                }
        # define material laws  - a tuple of dictionary
        self.materials = ({
                #material name
                "Name": "Mat1",
                #density
                "Density": 9.77e-4,
                #material type
                "Type": "Linear-Poroelastic",
                #"List of Elements": [(3,1),]
                # parameter of the law, each component of this list corresponds a line in the input file
                "Parameters": [(21,0.35,1,"PlaneStrain",),],
                # extra dof field
                "ExtraDof": ("DarcyLaw",(1.0e-6, 2.67e-5, 33.3, 1)),
                # support All or [(dim,phys), ...]
                "List of Elements": "All",
                },
                )
        self.initialBCs = {
                #"INITIAL CONDITIONS EXTRADOF": (0.,)
                }

    def getFEMConfigure(self):
        return self.FEM
    
    def getMMConfigure(self):
        return self.MM
    
    def getSolversConfigure(self):
        return self.solvers
        
    def getNeumannBCsConfigure(self):
        return self.NeumannBCs
    
    def getMaterialsConfigure(self):
        return self.materials 
    
    def getInitialBCsConfigure(self):
        return self.initialBCs
gmsh.initialize()

gmsh.open("box.msh")
model = gmsh.model()
#print(model.getDimension())
gmshModel = GmshModel.GmshModel(model)
config = myConfigure()

CreateInputFile.createInputFile(fileName="input.inp",model=gmshModel,configure=config)

gmsh.finalize()
