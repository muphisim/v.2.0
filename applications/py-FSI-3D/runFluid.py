import numpy as np
import pandas as pd
import os
import time 
import subprocess 
import matplotlib.pyplot as plt
import utils
from datetime import datetime

def writeGeoFile(boundFiles, geoFileName):
    nodeFiles = [f for f in boundFiles if "Node" in f]
    faceFiles = [f for f in boundFiles if "Face" in f]
    utils.createGeoFromBoundaryNodes(nodeFiles,faceFiles,geoFileName+"copy")
    inFile = open(geoFileName+"copy", 'r')
    outFile = open(geoFileName,"w")
    inTemplate = open("fluid.geo", 'r')
    for lines in inTemplate.readlines():
        outFile.write(lines)
        if lines.startswith("// Start solid"):
            for lines in inFile.readlines():
                if lines.startswith("SetFactory"):continue
                outFile.write(lines)
                if "flap = newv" in lines: break 
    outFile.close()

def extractFluidData(fluidFile, csvFileName, varNames: list, patchKeyName="flap"):
    if os.path.exists(csvFileName):df = pd.read_csv(csvFileName)
    else: df = pd.DataFrame({})
    f = open(fluidFile, 'r').readlines()
    for i, line in enumerate(f):
        if patchKeyName in line: ind = i
    try:        numCells = int(f[ind+4])
    except:        print(f[ind+4])
    if len(varNames) == 1:
        varList = []
        for i in range(numCells):
            try: varList.append(float(f[ind+6+i]))
            except: 
                print(f[ind+6+i])
                continue
        try: df[varNames[0]] = varList
        except: print(varList)
    elif len(varNames) == 3:
        varX, varY, varZ = [],[],[]
        for i in range(numCells):
            line = list(filter(lambda s: s.strip(), f[ind+6+i].split()))
            try:  
                varX.append(float(line[0][1:]))
                varY.append(float(line[1]))
                varZ.append(float(line[2][:-1]))
            except: 
                print(f[ind+6+i])
                continue
        try:
            df[varNames[0]] = varX
            df[varNames[1]] = varY
            df[varNames[2]] = varZ
        except: print(varX, varY, varZ)
    elif len(varNames) == 6:
        varXX, varXY, varXZ, varYY, varYZ, varZZ = [],[],[],[],[],[]
        for i in range(numCells):
            line = list(filter(lambda s: s.strip(), f[ind+6+i].split()))
            try:  
                varXX.append(float(line[0][1:]))
                varXY.append(float(line[1]))
                varXZ.append(float(line[2]))
                varYY.append(float(line[4]))
                varYZ.append(float(line[5]))
                varZZ.append(float(line[8][:-1]))
            except: 
                print(f[ind+6+i])
                continue
        try:
            df[varNames[0]],df[varNames[1]],df[varNames[2]],df[varNames[3]],df[varNames[4]],df[varNames[5]] = varXX, varXY, varXZ, varYY, varYZ, varZZ 
        except: print(varXX,varXY,varYY)
    else: print(varNames)
    df.drop(df.filter(regex='Unnamed').columns, axis=1, inplace=True) 
    df.to_csv(csvFileName)

def runFluid(caseDir, iteration, inBoundaryFiles, params, willDelete = False):
    """"
    Args:
        iteration (int) : iteration index
        inBoundaryFiles (list of strings) : boundary nodes provided by a list of in csv file names
        pressureFileName: distribution of pressure over the boundary is saved by this name
        parallel (boolean): True if the simulation is run in parallel
        endTime (int): maximum number of iterations in the simulation
        willDate(boolean): True if the running folder is deleted after finishing
        densityFluid (float): density of air to compute the force coefficients
    Returns:
        maxInt(float): the converged iteration 
    """
    boundFiles = [os.path.join(caseDir, ff) for ff in inBoundaryFiles]
    directory = os.getcwd()

    with open("timeLogFile.txt", "a") as f: f.write("\n starting to mesh fluid side :  "+str(datetime.now().strftime("%Y-%m-%d %H-%M-%S")))
    caseName = os.path.join(directory, f"{caseDir}/fluid-run-{iteration}")
    refCase = os.path.join(directory, f"{caseDir}/referenceFluid")
    os.system(f"rm -rf {caseName}")
    os.system(f"cp -rf {refCase} {caseName}")
    os.system(f"touch {caseName}/fluid-{iteration}.foam")
    geoFile = f"{caseName}/fluid.geo"
    mshFile = f"{caseName}/fluid.msh"
    #iteration = int(caseName.split("-")[-1])
    writeGeoFile(boundFiles, geoFile)
    subprocess.run(["gmsh","-3", f"{geoFile}", f"{mshFile}", "-v", "0","-nt", "4"])
    with open("timeLogFile.txt", "a") as f: f.write("\n Finished meshing, starting gmshToFoam : "+str(datetime.now().strftime("%Y-%m-%d %H-%M-%S")))
    os.system(f"gmshToFoam {mshFile} -case {caseName}  > {caseName}/runLog.txt;")
    with open("timeLogFile.txt", "a") as f: f.write("\n Finished gmshToFoam, starting polyMesh : "+str(datetime.now().strftime("%Y-%m-%d %H-%M-%S")))
    os.system(f"polyDualMesh 60 -overwrite -case {caseName}  >{caseName}/runLog.txt")
    with open("timeLogFile.txt", "a") as f: f.write("\n Finished polyMesh starting combineFaces : "+str(datetime.now().strftime("%Y-%m-%d %H-%M-%S")))
    os.system(f"checkMesh -case {caseName}  >{caseName}/checkMesh.txt")

    fluidProc = list(params["fluidProcessors"])
    totalProc = np.prod(fluidProc)
    with open("timeLogFile.txt", "a") as f: f.write(f"\n Finish checking mesh, starting fluid sim with {totalProc} processors: "+str(datetime.now().strftime("%Y-%m-%d %H-%M-%S")))
    #Write certain faces to empty
    filepath = os.path.join(caseName, "constant/polyMesh/boundary")
    boundary = open(filepath)
    boundary_list = boundary.readlines()
    ind1 = boundary_list.index("    flap\n")
    boundary_list[ind1+2] = "        type            wall;\n"
    boundary_list[ind1+3] = "        physicalType    wall;\n"
    boundary.close()
    with open(filepath, 'w') as f:
        f.writelines(boundary_list)

    # Write in the fluid processors decomposition
    filepath = os.path.join(caseName, "system/decomposeParDict")
    procs = open(filepath)
    procs_list = procs.readlines()
    for line in procs_list: 
        if "numberOfSubdomains  " in line: ind1 = procs_list.index(line)
        if "coeffs" in line: ind2 = procs_list.index(line)
    procs_list[ind1] = f"numberOfSubdomains  {totalProc};\n"
    procs_list[ind2+2] = f"        n            ({fluidProc[0]} {fluidProc[1]} {fluidProc[2]});\n"
    procs.close()
    with open(filepath, 'w') as f:
        f.writelines(procs_list)

    # Start running simulation
    if totalProc > 1:
        os.system(f"decomposePar -case {caseName} >{caseName}/log-decomposePar.txt")
        os.system(f"mpiexec -np {totalProc} simpleFoam -parallel -case {caseName} >{caseName}/runLog.txt")
        os.system(f"reconstructPar -case {caseName} >{caseName}/runLog.txt")
        os.system(f"rm -rf {caseName}/processor*")
    else:
        os.system(f"simpleFoam -case {caseName} >{caseName}/runLog.txt")
    with open("timeLogFile.txt", "a") as f: f.write("\n finish running fluid sim : "+str(datetime.now().strftime("%Y-%m-%d %H-%M-%S")))
    print(f"done running {caseName}")
    #time.sleep(0.1)
    maxInt = np.max([int(f) for f in os.listdir(caseName) if f.isdigit()])

    # Extract data from timepoint files
    fullCSV = os.path.join(f"{caseDir}/fluid-run-{iteration}","stress.csv")
    for file, var in zip(["C", "static(p)", "stressTensor"],
                        [["x", "y", "z"], ["p"], ["Sxx", "Sxy", "Sxz", "Syy", "Syz", "Szz"]]):
        extractFluidData(os.path.join(f"{caseDir}/fluid-run-{iteration}",str(maxInt), file), fullCSV, var)

    if willDelete:
        os.system(f"rm -rf {caseName}")

    return maxInt
