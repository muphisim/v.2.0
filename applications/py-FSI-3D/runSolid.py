import utils
import os
import CreateInputFile
import GmshModel
import gmsh
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

def nearestNeighbour(x, y, z, coords):
    """
    For point coords = (x0, y0, z0) obtain the closest point in the arrays x, y, z,
    and the position in the array of this closest point.
    """
    npoints = len(x)
    dist = np.zeros(npoints)
    for i in range(npoints):
        dist[i] =  (x[i] - coords[0])**2 + \
                   (y[i] - coords[1])**2 + \
                   (z[i] - coords[2])**2
    pos = dist.argmin()
    return pos, dist[pos]

def computeStress(gmshModel, pressureFileName, pressureBoundaries, boundaryFileNames):
    """
    Args:
        gmshModel: a model
        pressureFileName (string): path of the file storing pressure data
        pressureBoundaries (list of lists): boundaries where the pressure is applied [(dim, phys), ...]
        boundaryFileNames taken from previous solid iteration or initial Geo file.
    Return:
        Maps of facets and applied pressure, which is ready to include in the inp file
    """
    pData = pd.read_csv(pressureFileName)
    SXX = pData["Sxx"].values
    SXY = pData["Sxy"].values
    SXZ = pData["Sxz"].values
    SYY = pData["Syy"].values
    SYZ = pData["Syz"].values
    SZZ = pData["Szz"].values
    p = pData["p"].values
    pX = pData["x"].values
    pY = pData["y"].values
    pZ = pData["z"].values
    dData = pd.concat([pd.read_csv(f) for f in boundaryFileNames], axis=0)   
    origX = dData["origX"].values
    origY = dData["origY"].values
    origZ = dData["origZ"].values
    allEdges = {}
    for nn in pressureBoundaries:
        # pressureBoundaries is (1, 1), the FSI interface
        a = gmshModel.getGroupOfFacets(tuple(nn)) # returns [[1, [list of node points on interface 1]]]
        for kk in a:
            faceID = kk[0] # here is 1 for the FSI boundary index
            eleList = kk[1]
            for ele in eleList:
                face = GmshModel.GmshModel.getFacets(gmshModel.elementType,gmshModel.elements[ele-1])
                surfaceFace = face[faceID-1]
                surfaceNodes = list(set([node for facet in face for node in facet]))
                aveCentre = np.zeros(3)
                for node in surfaceFace: # face ID is 1 for this FSI so takes first entry
                    aveCentre += np.array(gmshModel.nodes[node-1])
                aveCentre *= (1./len(surfaceNodes))
                currLoc, _ = nearestNeighbour(origX, origY, origZ, aveCentre)
                interX, interY, interZ = dData.iloc[currLoc]["X"],dData.iloc[currLoc]["Y"],dData.iloc[currLoc]["Z"]
                pLoc, _ = nearestNeighbour(pX, pY, pZ, (interX, interY, interZ))
                # Take the stress tensor at this point and multiply by new solid normal
                T = np.array([[-p[pLoc] + SXX[pLoc], SXY[pLoc], SXZ[pLoc]],[SXY[pLoc], -p[pLoc] + SYY[pLoc], SYZ[pLoc]],[SXZ[pLoc], SYZ[pLoc], -p[pLoc] +SZZ[pLoc]]])#*SArea[pLoc]
                allEdges[(ele,faceID)] = [T[0][0], T[0][1], T[0][2], T[1][1], T[1][2], T[2][2]] # opposite convention in muphisim
    return allEdges
    

def runSolid(caseDir, iteration, initialBoundaryFileNames, currBoundaryFileNames, params):
    try:
        stressFileName = os.path.join(caseDir,f"fluid-run-{iteration}","stress.csv")
    except:
        print(f"File {stressFileName} does not exist")
        return None, None, np.nan

    runFolder = f"{caseDir}/solid-run-{iteration}"
    os.system(f"rm -rf {runFolder}")
    os.system(f"mkdir {runFolder}")
    outGeoFile = os.path.join(runFolder,"model.geo")
    outMshFile = os.path.join(runFolder,"solid.msh")
    boundFiles = [os.path.join(caseDir,ff) for ff in initialBoundaryFileNames]
    nodeFiles = [f for f in boundFiles if "Node" in f]
    faceFiles = [f for f in boundFiles if "Face" in f]
    with open("timeLogFile.txt", "a") as f: f.write("\n starting solid mesh : "+str(datetime.now().strftime("%Y-%m-%d %H-%M-%S")))
    
    if iteration==0:
        if os.path.exists(os.path.join(caseDir,  "referenceSolid","solid.msh")):
            print("trying to open meshFIle")
            prevMesh = os.path.join(caseDir,  "referenceSolid","solid.msh")
            os.system(f"cp {prevMesh} {outMshFile}")
            print("successfully opened meshFile")
        else:
            print("no reference solid mesh found")
            utils.createGeoFromBoundaryNodes(nodeFiles, faceFiles,outGeoFile)
            os.system(f"gmsh -{params['dimension']} {outGeoFile} -o {outMshFile} -v 0")
    else:
        prevMesh = os.path.join(caseDir, f"solid-run-{iteration-1}","solid.msh")
        os.system(f"cp {prevMesh} {outMshFile}")
    
    with open("timeLogFile.txt", "a") as f: f.write("\n finish meshing, start configuring solid : "+str(datetime.now().strftime("%Y-%m-%d %H-%M-%S")))
    gmsh.initialize()
    gmsh.open(outMshFile)
    model = gmsh.model()
    gmshModel = GmshModel.GmshModel(model)
    with open("timeLogFile.txt", "a") as f: f.write("\n starting to extract pressure: "+str(datetime.now().strftime("%Y-%m-%d %H-%M-%S")))
    stressBoundaries  = [[params["dimension"]-1, b] for b in params["fluidStressInterfaces"]]
    stressData = computeStress(gmshModel, stressFileName, stressBoundaries, boundaryFileNames = [os.path.join(caseDir,ff) for ff in currBoundaryFileNames if "Node" in ff])# This is very much a nearest neighbour mapping
    with open("timeLogFile.txt", "a") as f: f.write("\n start configuring: "+str(datetime.now().strftime("%Y-%m-%d %H-%M-%S")))
    
    materialProperties = [params["youngsModulus"],params["poissonRatio"],params["solidDensity"],params["solver"],params["solverNumSteps"]]
    config = utils.Configure(stressData, materialProperties)
    CreateInputFile.createInputFile(fileName=os.path.join(runFolder,"input.inp"),model=gmshModel,configure=config)
    gmsh.finalize()
    with open("timeLogFile.txt", "a") as f: f.write("\n finish configuring solid : "+str(datetime.now().strftime("%Y-%m-%d %H-%M-%S")))
    
    curDir = os.getcwd()
    print("ABOUT TO RUN MUPHISIM")
    nProcs = params["solidProcessors"]
    if nProcs ==1:os.system(f"MuPhiSim -input input.inp -inputDir {runFolder} -outputDir {runFolder} > {caseDir}/log-run.txt")
    else: os.system(f"mpirun -np {nProcs} MuPhiSim-mpi -input input.inp -inputDir {runFolder} -outputDir {runFolder} > {caseDir}/log-run.txt")
    print("DONE SIMULATION")
    os.chdir(curDir)
    with open("timeLogFile.txt", "a") as f: f.write("\n finish running solid sim : "+str(datetime.now().strftime("%Y-%m-%d %H-%M-%S")))
    # Find the last iteration number in the solid run
    allFiles = [d for d in os.listdir(runFolder) if d.startswith("MuPhiSimOutputs")]
    allNums = []
    for f in allFiles: allNums.append(float(f.split(".")[-2]))
    lastOutput = int(max(allNums))
    # Populate a dataframe with the unstressed position and displacement data from the output .vtk
    dispDF = {"x":[], "y":[], "z":[],"dx":[], "dy":[], "dz":[] }
    for outputFile in [f for f in os.listdir(runFolder) if f.endswith(f".{lastOutput}.vtk")]:
        with open(os.path.join(runFolder, outputFile), "r") as f: allLines = f.readlines()
        for i, line in enumerate(allLines):
            if line.startswith("POINTS "): pointIdx = i
            if line.startswith("VECTORS Displacement float"): dispIdx = i
        for i, line in enumerate(allLines):
            if i < pointIdx+1: continue
            try:
                dispDF["x"].append(float(line.split(" ")[0]))
                dispDF["y"].append(float(line.split(" ")[1]))
                dispDF["z"].append(float(line.split(" ")[2]))
            except: break
        for i, line in enumerate(allLines):
            if i < dispIdx+1: continue
            try:
                dispDF["dx"].append(float(line.split(" ")[0]))
                dispDF["dy"].append(float(line.split(" ")[1]))
                dispDF["dz"].append(float(line.split(" ")[2]))
            except: break
    # save this data for reference
    dispDF = pd.DataFrame(dispDF) 
    dispDF.to_csv(os.path.join(runFolder, "extractedData.csv"), index=False)

    maxDisp, minDistShunt = -1, 100
    nodeFiles = [os.path.join(caseDir,ff) for ff in currBoundaryFileNames if "Node" in ff]
    for ff in nodeFiles:
        newLines, ids = [],[]
        lastDf = pd.read_csv(ff)
        outfile = os.path.join(runFolder,ff.split("/")[-1])
        for index, row in lastDf.iterrows():  # row is previous position, needs to match up to new one
            centre = (row["origX"], row["origY"], row["origZ"]) 
            # unstressed position of node for previous data should match to current 
            pos, dist = nearestNeighbour(dispDF["x"], dispDF["y"], dispDF["z"], centre)
            if dist < 1e-4: 
                # calculate current position for this iteration
                newPos = np.array([dispDF.at[pos, "x"]+dispDF.at[pos, "dx"],dispDF.at[pos, "y"]+dispDF.at[pos, "dy"],dispDF.at[pos, "z"]+dispDF.at[pos, "dz"]])
                unstressedPos = np.array([dispDF.at[pos, "x"],dispDF.at[pos, "y"],dispDF.at[pos, "z"]])
                # position from previous iteration
                lastPos = np.array([row["X"], row["Y"], row["Z"]])
                dX = params["relaxationModulus"]*(newPos-lastPos)
                newPos = lastPos + dX
                maxDisp =np.max([maxDisp, np.linalg.norm(newPos-unstressedPos)])
                ids.append(int(row["nodeID"]))
                newLines.append(list(newPos) + list(dX) + [row["origX"], row["origY"], row["origZ"]])
        p = pd.DataFrame()
        p["nodeID"] = ids
        p[["X","Y","Z","dX","dY","dZ","origX","origY","origZ"]] = newLines
        p.to_csv(outfile, index=False)  
        mD = utils.minDistToShunt(outfile)
        if minDistShunt > mD: minDistShunt = mD 
    for ff in faceFiles:   
        outFile = os.path.join(runFolder,ff.split("/")[-1])
        os.system(f"cp {ff} {outFile}")
    return maxDisp,  minDistShunt
