import gmsh
import GmshModel
import os
import utils
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

os.system("rm -rf initialGeo")
os.system("mkdir initialGeo")

utils.createGeo("0012",0.5,"initialGeo/model.geo")
os.system("gmsh -2 -order 2 initialGeo/model.geo -o initialGeo/model.msh")

gmsh.initialize()
gmsh.open("initialGeo/model.msh")
model = gmsh.model()
gmshModel = GmshModel.GmshModel(model)

plt.figure()

nodeGrps = [(1,3), (1,4), (1,5)]
allFiles = []
allSaveNodes = []
for oneGr in nodeGrps:
    saveNodes = utils.getBoundaryNodes(gmshModel,[oneGr])
    allSaveNodes.append(saveNodes)
    nn = [n-1 for n in saveNodes]
    boundaryNodes= gmshModel.nodes[nn]
    plt.plot(boundaryNodes[:,0],boundaryNodes[:,1],".-",label=f"boundaryNodes{nodeGrps.index(oneGr)}")
    saveFileName = f"initialGeo/boundaryNodes{nodeGrps.index(oneGr)}.csv"
    allFiles.append(saveFileName)
    p = pd.DataFrame()
    p["nodeID"] = np.array(saveNodes)
    p[["X","Y","Z"]] = boundaryNodes
    p.to_csv(saveFileName,index=False)

plt.xlabel("x")
plt.ylabel("y")
plt.axis('equal')
plt.legend()

plt.savefig("initialGeo/initialProfile.png",dpi=1000)

utils.createGeoFromBoundaryNodes(allFiles,"initialGeo/myGeo.geo")


