from Aerofoil import *
import pandas as pd
import numpy as np
import gmsh
import GmshModel
import Configure
import os

def createGeo(nacaName, xs, geoFileName):
    """
    Args:
        nacaName:   Name of Naca 4 digits
        xcMin:      min xc
        xcMax:      max of xc 
        geoFileName:    output file name with gmsh
    Returns:
        None
    """
    upper, lower = MPXX_family(nacaName,normalised = False)
    xu = np.flip(upper[:,0])
    yu = np.flip(upper[:,1])
    xl = lower[1:-1,0]
    yl = lower[1:-1,1]
    numU= len(xu)
    numL = len(xl)
    locPu = 0
    locPl = 0
    for i in range(numU):
        if xu[i] < xs:
            locPu = i+1
            break
    for i in range(numL):
        if xl[i] > xs:
             locPl = i+numU
             break
   
    file = open(geoFileName,"w")
    file.write("General.ExpertMode = 1;\n")
    file.write("lsca=0.005;\n")
    for ii in range(numU):
        file.write(f"Point({ii+1})={{{xu[ii]},{yu[ii]},0,lsca}};\n")
    for ii in range(numL):
        file.write(f"Point({ii+1+numU})={{{xl[ii]},{yl[ii]},0,lsca}};\n")
        
    file.write(f"ls1 = newl; Spline(ls1)={{1:{locPu}}};\n")
    file.write(f"ls2 = newl; Spline(ls2)={{{locPu}:{locPl}}};\n")
    file.write(f"ls3 = newl; Spline(ls3)={{{locPl}:{numU+numL},1}};\n")
    file.write(f"ls4 = newl; Line(ls4)={{{locPu},{locPl}}};\n")
    file.write("loop1=newl;Line Loop(loop1)={ls1,ls3,ls4};\n")
    file.write("loop2=newl;Line Loop(loop2)={-ls4,ls2};\n")
    file.write("s1 = news; Plane Surface(s1)={loop1};\n")
    file.write("s2 = news; Plane Surface(s2)={loop2};\n")
    file.write("Physical Surface(1) = {s1};\n")
    file.write("Physical Surface(2) = {s2};\n")
    file.write("Physical Curve(\"leading\",3) = {ls2};\n")
    file.write("Physical Curve(\"trailingUp\",4) = {ls1};\n")
    file.write("Physical Curve(\"trailingLow\",5) = {ls3};\n")
    file.close()
    
def getBoundaryNodes(gmshModel, nodeGrps):
    """
    Args:
        gmshModel:   gmsh model
        grp: group index by vector of tuple of two elements, eg [(dim, phys),]
    Returns:
        NodeList

    """
    elements = gmshModel.elements
    edges= []
    for nn in nodeGrps:
        a = gmshModel.getGroupOfFacets(nn)
        for kk in a:
            faceID = kk[0]
            eleList = kk[1]
            for ele in eleList:
                line = GmshModel.GmshModel.getEdges(gmshModel.elementType,elements[ele-1])
                edges.append(line[faceID-1])
    nodeList = [edges[0][0]]
    if len(edges[0]) > 2:
        for k in edges[0][2:]:
            nodeList.append(k)
    nodeList.append(edges[0][1])
    newedges = list(edges[1:])
    ite=0
    while True:
        ite = ite +1
        found = None
        for edId in range(len(newedges)):
            ed = newedges[edId]
            if ed[0] == nodeList[-1]:
                nodeList.append(ed[1])
                found = edId
                break
            elif ed[1] == nodeList[-1]:
                nodeList.append(ed[0])
                found = edId
                break
            elif ed[0] == nodeList[0]:
                nodeList = [ed[1]]+ nodeList
                found = edId
                break
            elif ed[1] == nodeList[0]:
                nodeList = [ed[0]]+ nodeList
                found = edId
                break
        if found is not None:
            newedges.pop(found)
        if len(newedges) == 0:
            break
        
        if ite > 100000:
            raise RuntimeError("maximal number of iterations is reached!")
        
    if nodeList[0] == nodeList[-1]:
        return nodeList[:-1]
    else:
        return nodeList

def getDeformedNodesRelaxation(allNodes, saveNodes, allDisp, omega=1., allPosPrev=None, saveFileName=None):
    newNodes = np.array(allNodes)
    maxChange = 0.
    if allPosPrev is None:
        for node in saveNodes:
            dx = allDisp[f"Node{node}Comp0"]
            dy = allDisp[f"Node{node}Comp1"]
            newNodes[node-1][0] += omega*dx
            newNodes[node-1][1] += omega*dy
            if maxChange < np.abs(omega*dx):
                maxChange = np.abs(omega*dx)
            if maxChange < np.abs(omega*dy):
                maxChange = np.abs(omega*dy)
    else:
        for node in saveNodes:
            dx = allDisp[f"Node{node}Comp0"]
            dy = allDisp[f"Node{node}Comp1"]
            xPrev = allPosPrev[f"Node{node}Comp0"]
            yPrev = allPosPrev[f"Node{node}Comp1"]
            posX = newNodes[node-1][0]
            posY = newNodes[node-1][1]
            newNodes[node-1][0] = posX + (xPrev-posX)+ omega*(dx-(xPrev-posX))
            newNodes[node-1][1] = posY + (yPrev-posY)+ omega*(dy-(yPrev-posY))
            
            if maxChange < np.abs(newNodes[node-1][0]-posX):
                maxChange = np.abs(newNodes[node-1][0]-posX)
            if maxChange < np.abs(newNodes[node-1][1]- posY ):
                maxChange = np.abs(newNodes[node-1][1]- posY )
            
    nn = [n-1 for n in saveNodes]
    bdNodes= newNodes[nn]
    if saveFileName is not None:
        p = pd.DataFrame()
        p["nodeID"] = np.array(saveNodes)
        p[["X","Y","Z"]] = bdNodes
        p.to_csv(saveFileName,index=False)
    return maxChange

def createProfile(boundaryNodeFiles):
    """
    create profile from boundary nodes 
    """
    nodeData = []
    nodeIds = {}
    nodePerGroups = {}
    for ff in boundaryNodeFiles:
        p = pd.read_csv(ff)   
        nodeData.append(p.values.tolist())
        nodePerGroups[boundaryNodeFiles.index(ff)] = []
        for n in nodeData[-1]:
            nodePerGroups[boundaryNodeFiles.index(ff)].append(int(n[0]))
            if n[0] not in nodeIds.keys():
                nodeIds[int(n[0])] = n[1:]

    interSecNode = list(set(nodePerGroups[1])&set(nodePerGroups[2]))
    pos1 = nodePerGroups[1].index(interSecNode[0])
    pos2 = nodePerGroups[2].index(interSecNode[0])
    if pos1 > 0:
        nodePerGroups[1].reverse()
    if pos2 ==0:
        nodePerGroups[2].reverse()
     
    if nodePerGroups[0][0] != nodePerGroups[1][-1]:
        nodePerGroups[0].reverse()
        
    return np.array([nodeIds[i] for i in nodePerGroups[1]]+
                    [nodeIds[i] for i in nodePerGroups[0][1:]]+
                    [nodeIds[i] for i in nodePerGroups[2][1:-1]])
        
        
    
def createGeoFromBoundaryNodes(boundaryNodeFiles, geoFileName):
    """
    Args:
        boundaryNodeFiles:   Vector of csv files
    Returns:
        None
    """
    nodeData = []
    nodeIds = {}
    nodePerGroups = {}
    for ff in boundaryNodeFiles:
        p = pd.read_csv(ff)   
        nodeData.append(p.values.tolist())
        nodePerGroups[boundaryNodeFiles.index(ff)] = []
        for n in nodeData[-1]:
            nodePerGroups[boundaryNodeFiles.index(ff)].append(int(n[0]))
            if n[0] not in nodeIds.keys():
                ss = len(nodeIds)+1
                nodeIds[int(n[0])] = [ss, n[1:]]
    #print(nodeIds)
    #print(nodePerGroups)
    
    interSecNode = list(set(nodePerGroups[1])&set(nodePerGroups[2]))
    pos1 = nodePerGroups[1].index(interSecNode[0])
    pos2 = nodePerGroups[2].index(interSecNode[0])
    if pos1 > 0:
        nodePerGroups[1].reverse()
    if pos2 >0:
        nodePerGroups[2].reverse()
     
    if nodePerGroups[0][0] != nodePerGroups[1][-1]:
        nodePerGroups[0].reverse()
     
 
    file = open(geoFileName,"w")
    file.write("General.ExpertMode = 1;\n")
    file.write("lsca=0.0075;\n")
    for nn, vv in nodeIds.items():
        nodeIndex = vv[0]
        x= vv[1][0]
        y= vv[1][1]
        
        file.write(f"Point({nodeIndex})={{{x},{y},0,lsca}};\n")
        
    for lineIndex, nodes in nodePerGroups.items():
        for ii in range(1,len(nodes)):
            first = nodeIds[nodes[ii-1]][0]
            second = nodeIds[nodes[ii]][0]
            file.write(f"ls{lineIndex}[{ii-1}] = newl; Line(ls{lineIndex}[{ii-1}])={{{first},{second}}};\n")
      
    file.write(f"li = newl; Line(li)={{{nodeIds[nodePerGroups[0][0]][0]},{nodeIds[nodePerGroups[0][-1]][0]}}};\n")

    file.write("loop1=newl;Line Loop(loop1)={ls0[],-li};\n")
    file.write("s1 = news; Plane Surface(s1)={loop1};\n")
    file.write("Physical Surface(\"Part1\",1) = {s1};\n")
    
    file.write("loop2=newl;Line Loop(loop2)={ls1[],li,-ls2[]};\n")
    file.write("s2 = news; Plane Surface(s2)={loop2};\n")
    file.write("Physical Surface(\"Part2\",2) = {s2};\n")
 
    file.write("Physical Curve(\"leading\",3) = {ls0[]};\n")
    file.write("Physical Curve(\"trailingUp\",4) = {ls1[]};\n")
    file.write("Physical Curve(\"trailingLow\",5) = {ls2[]};\n")
    file.close()
    return (1,3), (1,4), (1,5)
    
class myConfigure (Configure.Configure):
    def __init__(self, saveNodes, pressureData, materialProperties):
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
                #scale factor to estimate time sep for solver
                "Number Of Steps": 1,
                #output specification, each component of the list corresponds a line in the oxfemm input file
                "Outputs": ((1000,), ("All" , "max"),),
                # start and end time
                "Time": (0.,1.),
                "DispBC":(
                            # each bc contains (type, location(dim, physical) list of prescried components and prescribed values
                            ("DISPLACEMENT RAMP",(2,1), ((0,0.),(1,0.))), 
                         ),
                "Absolute Tolerance": 1e-6,
                "Relative Tolerance": 1e-6,
                "PETSC_SOLVER":"lu",
                "Extraction":[
                                
                               ("NODE",f"Group3","Unknown","Rough",(1,3),0),
                               ("NODE",f"Group4","Unknown","Rough",(1,4),0),
                               ("NODE",f"Group5","Unknown","Rough",(1,5),0),
                               ("NODE",f"Group3","Unknown","Rough",(1,3),1),
                               ("NODE",f"Group4","Unknown","Rough",(1,4),1),
                               ("NODE",f"Group5","Unknown","Rough",(1,5),1),
                               ("NODE",f"Group3","ExternalForce","Rough",(1,3),0),
                               ("NODE",f"Group4","ExternalForce","Rough",(1,4),0),
                               ("NODE",f"Group5","ExternalForce","Rough",(1,5),0),
                               ("NODE",f"Group3","ExternalForce","Rough",(1,3),1),
                               ("NODE",f"Group4","ExternalForce","Rough",(1,4),1),
                               ("NODE",f"Group5","ExternalForce","Rough",(1,5),1),
                            ]
                },
                )
        # define all Neumann BCs
        self.NeumannBCs = {
                "PressureBC": (
                                # type, location, time interval, value
                                #("PRESSURE RAMP", (1,4), (0.,10.), 1e3),
                                #("PRESSURE RAMP", (2,52), (0.5,1.), -0.1),
                              ),
                "ElementPressureBC": [
                                        ("PRESSURE RAMP", (edge[0],edge[1]), (0.,10.), p) 
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
                "Density": 1.,
                #material type
                "Type": "Linear-Elastic",
                # parameter of the law, each component of this list corresponds a line in the input file
                "Parameters": [(materialProperties[0],materialProperties[1],"PlaneStrain"),],
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
