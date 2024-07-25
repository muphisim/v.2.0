import utils
import os
import CreateInputFile
import GmshModel
import gmsh
import pandas as pd
import numpy as np


def getPressureAtPointFromData(x, y, z, p, coords):
    """
    get pressure at point coords = (x0, y0, z0) from data specified by (x, y, z, p)
    """
    npoints = len(p)
    dist = np.zeros(npoints)
    for i in range(npoints):
        dist[i] =  (x[i] - coords[0])**2 + \
                   (y[i] - coords[1])**2 + \
                   (z[i] - coords[2])**2
    pos = dist.argmin()
    return p[pos]


def computePressure(gmshModel, pressureFileName, pressureBoundaries):
    """
    Args:
        gmshModel: a model
        pressureFileName (string): path of the file storing pressure data
        pressureBoundaries (list of lists): boundaries where the pressure is applied [(dim, phys), ...]
    Return:
        Maps of facets and applied pressure, which is ready to include in the inp file
    """
    print("start computePressure")
    pData = pd.read_csv(pressureFileName)
    x = pData.values[:,0]
    y = pData.values[:,1]
    z = pData.values[:,2]
    p = pData.values[:,4]
    allEdges = {}
    #allNodes = []
    for nn in pressureBoundaries:
        a = gmshModel.getGroupOfFacets(nn)
        for kk in a:
            faceID = kk[0]
            eleList = kk[1]
            for ele in eleList:
                line = GmshModel.GmshModel.getEdges(gmshModel.elementType,gmshModel.elements[ele-1])                
                aveCentre = np.zeros(3)
                for node in line[faceID-1]:
                    aveCentre += np.array(gmshModel.nodes[node-1])
                aveCentre *= (1./len(line[faceID-1]))
                ploc = getPressureAtPointFromData(x, y, z, p, aveCentre)
                #print(aveCentre, ploc)
                allEdges[(ele,faceID)] = -ploc # opposite convention in muphisim
                #allNodes.append(aveCentre.tolist()+[ploc])
    
    
    print("done computePressure")
    #return np.array(allNodes), pData.values
    return allEdges
    

def runSolid(caseDir, iteration, inBoundaryFileNames, pressureFileName, pressureBoundaries, materialProperties, relaxFactor=1.):
    
    if not(os.path.exists(pressureFileName)):
        print(f"File {pressureFileName} does not exist")
        return None, None, np.nan

    runFolder = f"{caseDir}/solid-run-{iteration}"
    os.system(f"rm -rf {runFolder}")
    os.system(f"mkdir {runFolder}")
    outGeoFile = os.path.join(runFolder,"model.geo")
    outMshFile = os.path.join(runFolder,"model.msh")
    boundFiles = [os.path.join(caseDir,ff) for ff in inBoundaryFileNames]
    nodeGrps = utils.createGeoFromBoundaryNodes(boundFiles,outGeoFile)
    os.system(f"gmsh -2 -order 2 {outGeoFile} -o {outMshFile}")
 
    gmsh.initialize()
    gmsh.open(outMshFile)
    model = gmsh.model()
    gmshModel = GmshModel.GmshModel(model)
    saveNodes = set()
    for nn in nodeGrps:
        nodes = gmshModel.getGroupOfNodes(nn)
        for ni in nodes:
            saveNodes.add(ni)
    saveNodes = list(saveNodes)
    pressureData = computePressure(gmshModel, pressureFileName,pressureBoundaries)
        
    config = utils.myConfigure(saveNodes, pressureData, materialProperties)
    CreateInputFile.createInputFile(fileName=os.path.join(runFolder,"input.inp"),model=gmshModel,configure=config)
    gmsh.finalize()
    
    curDir = os.getcwd()
    os.system(f"MuPhiSim -input input.inp -inputDir {runFolder} -outputDir {runFolder}")
    print("DONE SIMULATION")
    os.chdir(curDir)
    
    dataFiles = ["Group3_Unknown0_Rough.csv","Group3_Unknown1_Rough.csv",
            "Group4_Unknown0_Rough.csv","Group4_Unknown1_Rough.csv",
            "Group5_Unknown0_Rough.csv","Group5_Unknown1_Rough.csv",]
    allDisp = {}
    for ff in dataFiles:
        if not(os.path.exists(os.path.join(runFolder,ff))):
            print(f"File {os.path.join(runFolder,ff)} does not exist")
            return None, None, np.nan, np.nan
        data = pd.read_csv(os.path.join(runFolder,ff))
        for col in data.columns[1:]:
            allDisp[col] = data[col].values[-1]
                
    allPosPrev = None
    if iteration >0:
        allPosPrev={}
        inBoundaryFiles = [f"{caseDir}/solid-run-{iteration-1}/boundaryNodes0.csv",
                     f"{caseDir}/solid-run-{iteration-1}/boundaryNodes1.csv",
                     f"{caseDir}/solid-run-{iteration-1}/boundaryNodes2.csv"]
                     
        for ff in inBoundaryFiles:
            if not(os.path.exists(os.path.join(ff))):
                print(f"File {ff} does not exist")
                return None, None, np.nan, np.nan
            data = pd.read_csv(ff)
            for row in data.values:
                allPosPrev[f"Node{int(row[0])}Comp0"] = row[1]
                allPosPrev[f"Node{int(row[0])}Comp1"] = row[2]
                
            
    forceFiles = ["Group3_ExternalForce0_Rough.csv","Group3_ExternalForce1_Rough.csv",
            "Group4_ExternalForce0_Rough.csv","Group4_ExternalForce1_Rough.csv",
            "Group5_ExternalForce0_Rough.csv","Group5_ExternalForce1_Rough.csv",]
    
    maxForce = 0
    for ff in forceFiles:
        if not(os.path.exists(os.path.join(runFolder,ff))):
            print(f"File {os.path.join(runFolder,ff)} does not exist")
            return None, None, np.nan, np.nan
        data = pd.read_csv(os.path.join(runFolder,ff))
        for col in data.columns[1:]:
            if np.abs(data[col].values[-1]) > maxForce:
                maxForce = np.abs(allDisp[col])
    
    maxDisp = 0
    for ii in range(len(nodeGrps)):
        saveNodes = utils.getBoundaryNodes(gmshModel,[nodeGrps[ii]])
        fileSplit = inBoundaryFileNames[ii].split("/")
        outFile = os.path.join(runFolder,fileSplit[-1])
        maxChange = utils.getDeformedNodesRelaxation(gmshModel.nodes,saveNodes,allDisp,relaxFactor, allPosPrev ,outFile)
        print(maxChange)
        if maxDisp < maxChange:
            maxDisp = maxChange
   
    return gmshModel, pressureData, maxDisp, maxForce
   
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    plt.close("all")
   
    CpFile = "fluid-run-1000/Cp.csv"
    pressure = pd.read_csv(CpFile)
    pressure = pressure.values
    iteration=0
    inBoundaryFiles = ["initialGeo/boundaryNodes0.csv",
                 "initialGeo/boundaryNodes1.csv",
                 "initialGeo/boundaryNodes2.csv"]
    gmshModel, pressureData, maxDisp, maxForce = runSolid(".",iteration,inBoundaryFiles,CpFile,[(1,4),(1,5)],[5e4,0.3])
    nodeData = []
    for kk, p in pressureData.items():
        ele = kk[0]
        edgeId = kk[1]
        line = GmshModel.GmshModel.getEdges(gmshModel.elementType,gmshModel.elements[ele-1]) 
        aveCentre = np.zeros(3)
        for node in line[edgeId-1]:
            aveCentre += np.array(gmshModel.nodes[node-1])
        aveCentre *= (1./len(line[edgeId-1]))
        nodeData.append(aveCentre.tolist()+[p])
    pressureData = np.array(nodeData)
    
    plt.figure()
    plt.plot(pressureData[:,0],-pressureData[:,3],".",label="approx")
    plt.plot(pressure[:,0],pressure[:,4],".", label="correct")
    plt.legend()
    plt.show()
