import gmsh
import GmshModel
import os
import utils
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json
from datetime import datetime

params = json.loads(open("config.json", "r").read())
nodeGrps = [tuple([params["dimension"]-1, b]) for b in params["deformingInterfaces"]] # these tuples are dimension, index of the deformeable boundary??
    
def getBoundaryNodes(gmshModel, nodeGrps):
    """
    Args:
        gmshModel:   gmsh model
        grp: group index by vector of tuple of two elements, eg [(dim, phys),]
    Returns:
        NodeList

    """
    elements = gmshModel.elements
    faces= []
    for nn in nodeGrps:
        a = gmshModel.getGroupOfFacets(tuple(nn))
        for kk in a:
            faceID = kk[0]
            eleList = kk[1]
            for ele in eleList:
                face = GmshModel.GmshModel.getFacets(gmshModel.elementType,elements[ele-1])
                faces.append( face[faceID-1])
    faceNodes = list(set([node for facet in faces for node in facet]))
    return np.array(faces), faceNodes       
 
## Insert updated mesh size
dim = params["dimension"]
for meshSize, file in zip([params["fluidMeshSize"],params["solidMeshSize"]], ["fluid.geo", "solid.geo"]):#[j for j in os.listdir(".") if '.geo' in j]:
    with open(file, "r") as f: allLines = f.readlines()
    for i, line in enumerate(allLines):
        if "meshSize" in line: 
            allLines[i] = f"meshSize = {meshSize};\n"
            break
    with open(file, "w") as f: f.writelines(allLines)

os.system("rm -rf initialGeo")
os.system("mkdir initialGeo")
with open("timeLogFile.txt", "a") as f: f.write("Starting initial config :  "+str(datetime.now().strftime("%Y-%m-%d %H-%M-%S")))
os.system(f"gmsh -{params['dimension']} solid.geo -o initialGeo/solid.msh -v 0") # The option -order 2 causes problems

gmsh.initialize()
gmsh.open("initialGeo/solid.msh")
model = gmsh.model()
gmshModel = GmshModel.GmshModel(model)

allFiles, faceFiles = [], []
allSaveNodes, allSaveFaces = [], []
for oneGr in nodeGrps:
    faces, saveNodes = getBoundaryNodes(gmshModel,[oneGr])
    allSaveNodes.append(saveNodes)
    allSaveFaces.append(faces)
    nn = [n-1 for n in saveNodes]
    boundaryNodes= gmshModel.nodes[nn]
    saveFileName = f"initialGeo/boundaryNodes{nodeGrps.index(oneGr)}.csv"
    allFiles.append(saveFileName)
    p = pd.DataFrame()
    p["nodeID"] = np.array(saveNodes)
    p[["X","Y","Z"]] = boundaryNodes
    p[["origX","origY","origZ"]] = boundaryNodes
    p.to_csv(saveFileName,index=False)

    saveFileName = f"initialGeo/boundaryFaces{nodeGrps.index(oneGr)}.csv"
    faceFiles.append(saveFileName)
    p = pd.DataFrame()
    p["faceID"] = np.array([i for i in range(len(faces))])
    p[["node1","node2","node3"]] = faces
    p.to_csv(saveFileName,index=False)
utils.createGeoFromBoundaryNodes(allFiles,faceFiles,"initialGeo/solid.geo")
with open("timeLogFile.txt", "a") as f: f.write("\n finishing initial config :  "+str(datetime.now().strftime("%Y-%m-%d %H-%M-%S")))


