
import pandas as pd
import numpy as np
import Configure
import json
from itertools import product


params = json.loads(open("config.json", "r").read())

def createGeoFromBoundaryNodes(boundaryNodeFiles, boundaryFaceFiles, geoFileName):
    """
    Args:
        boundaryNodeFiles:   Vector of csv files
        boundaryFaceFiles:   Vector of csv files
    Returns:
        None
    """    
    params = json.loads(open("config.json", "r").read())
    file = open(geoFileName,"w")
    # Any precursor lines 
    file.write(f'SetFactory("OpenCASCADE");General.ExpertMode = 1;meshSize={params["solidMeshSize"]};\n')
    Lz, flapWidth = -1, -1

    nodeIds = []
    for ff in boundaryNodeFiles:
        p = pd.read_csv(ff)  
        newNode = p.values.tolist()
        for n in newNode:
            if n[0] not in nodeIds:
                nodeIds.append(n[0])
                if n[3] > Lz:Lz = n[3]
                if n[3]==0 and n[2] < 1e-4 and np.abs(n[1])>flapWidth:flapWidth = np.abs(n[1])
                file.write(f"Point({int(n[0])})={{{n[1]},{n[2]},{n[3]},meshSize}};\n")
    
    for ff in boundaryFaceFiles:
        p = pd.read_csv(ff)
        surfaceID = ff.split(".")[0][-1]
        for _, row in p.iterrows():
            faceID, pt1, pt2, pt3 = int(row["faceID"]), int(row["node1"]), int(row["node2"]),int(row["node3"])
            file.write(f"l = newl; Line(l)={{{pt1},{pt2}}};l = newl; Line(l)={{{pt2},{pt3}}};l = newl; Line(l)={{{pt3},{pt1}}};\n")
            file.write('ll = newll; Line Loop(ll)={l-2, l-1, l};')
            file.write(f"ls{surfaceID}[{int(faceID)}]=news; Plane Surface(ls{surfaceID}[{int(faceID)}])={{ll}};\n")

    file.write(f"bottom = news; Rectangle(bottom) = {{-{flapWidth}, 0, 0,{2*flapWidth}, {Lz}}};\n")
    file.write("Rotate{{1, 0, 0}, {0, 0, 0}, Pi/2}{Surface{bottom};}\n")
    file.write("Coherence;loop1=newsl;Surface Loop(loop1)={ls0[],ls1[], bottom};\n")
    file.write("Coherence;flap = newv; Volume(flap) = {loop1};\n")# Recombine Surface{s1};\n")
    file.write("Physical Surface(\"moving\", 1) = {ls0[]};\n") 
    file.write("Physical Surface(\"slip\", 2) = {ls1[]};\n") 
    file.write("Physical Surface(\"fixed\", 3) = {bottom};\n")
    file.write("Physical Volume(\"solid\", 1) = {flap};\n")
    file.close()
    return

def dist(pts):
    x1, y1, z1, x2, y2, z2 = pts[0][0], pts[0][1], pts[0][2], pts[1][0], pts[1][1], pts[1][2]
    return  np.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
def closeToCentre(pts):
    ## Return true if any deformed flap points have intersected with the catheter
    x1, y1, z1, xc, yc, r = pts[0][0], pts[0][1], pts[0][2], pts[1][0], pts[1][1], pts[1][2]
    if dist([[x1, y1, 0], [xc, yc, 0]]) < r: 
        return True
    else: return False

import matplotlib.pyplot as plt
def minDistToShunt(boundaryNodeFile):
    xc, yc, L, r, n = 2, 2.5, 30, 1.25, 50
    shuntPoints = {"x": [xc+r*np.cos(np.pi*i/n + np.pi/2) for i in range(n)],"y": [yc+r*np.sin(np.pi*i/n + np.pi/2) for i in range(n)],
                   "z": [L/n*i for i in range(n)]}
    shuntPoints = product([[x, y] for x, y in zip(shuntPoints["x"],shuntPoints["y"])], list(shuntPoints["z"]))
    shuntPoints = [(s[0][0], s[0][1], s[1]) for s in list(shuntPoints)]
    p = pd.read_csv(boundaryNodeFile)
    newX, newY, newZ = p["X"].values, p["Y"].values, p["Z"].values
    p = list(zip(newX, newY, newZ))
    if closeToCentre(max(product(p, [(xc, yc, r)]), key=closeToCentre)): return 0
    else:
        point_pairs = product(shuntPoints, p)
        minDisp = dist(min(point_pairs, key=dist))
        return minDisp
    
class Configure (Configure.Configure):
    def __init__(self, pressureData, materialProperties):
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
                "Solver": materialProperties[3][1:-1],
                #scale factor to estimate time sep for solver
                "Scale Factor": 1,
                "Number Of Steps": materialProperties[4],
                #output specification, each component of the list corresponds a line in the oxfemm input file
                "Outputs": ((10,), ("All" , "max"),),
                # start and end time
                "Time": (0.,100.),
                "DispBC":(
                            # each bc contains (type, location(dim, physical) list of prescried components and prescribed values
                            ("DISPLACEMENT RAMP",(2, 3), ((0,0.),(1,0.),(2,0.))), 
                            ("DISPLACEMENT RAMP",(2, 2), ((2,0.),)), 
                         ),
                "Absolute Tolerance": 1e-6,
                "Relative Tolerance": 1e-6,
                #"PETSC_SOLVER":"lu",
                #"Extraction":[("NODE",f"Group{i}","Unknown","Rough",(2,i),j) for i, j in zip(3*[1,2],len([1,2])*[0,1,2])]
                                #[("NODE",f"Group{i}","ExternalForce","Rough",(2,i),j) for i, j in  zip(3*[1,2],len([1,2])*[0,1,2])]
                },
                )

        # define all Neumann BCs
        self.NeumannBCs = {
                "ElementPressureBC": [
                                    #("PRESSURE INST", (edge[0],edge[1]), (0.,1.), p) 
                                    #("PRESSURE RAMP", (edge[0],edge[1]), (0.,1.), p[0]) 
                                     #   for edge, p in pressureData.items()
                                    ],
                "ElementStressBC": [
                                        ("STRESS RAMP", (edge[0],edge[1]), (0.,100.), (p[0], p[1],p[2], p[3], p[4], p[5])) 
                                        for edge, p in pressureData.items()
                                    ],
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
                "Density": materialProperties[2],
                #material type
                # parameter of the law, each component of this list corresponds a line in the input file
                "Type": "HyperElastic, Neo-Hookean",
                "Parameters": [(materialProperties[0],materialProperties[1],),],
                # support All or [(dim,phys), ...]
                #"ExtraDof": ("Temperature",(1.e-10, 0, 5, 0, 0, 0, 0)),
                "List of Elements": "All"
                #"List of Elements": [(3,1),]
                },
                )
        self.initialBCs = {
                #"INITIAL CONDITIONS EXTRADOF": (273.,)
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

