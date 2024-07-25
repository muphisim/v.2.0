import numpy as np
import pandas as pd
import os
from Aerofoil import *
import run2DCase  
import utils

def getCp(folderName,step,fileName, Uinf, densityAir, saveFile):
    """
    Args:
        folderName (string): postprocessing folder
        step (int) : step at which the solution converges
        fileName (string): the name for Cp
        Uinf (float) : farfield flow velocity magnitude
        densityAir (float) : density of air
        saveFile (string): the full-path of the file that pressure is stored
    Returns:
        x, y, z, Cp, p: coordinates, pressure coefficient and pressure (difference with respect to the air pressure) 
    """
    fname = os.path.join(folderName,str(step),fileName)
    value = np.genfromtxt(fname,comments="#")
    x = value[:,0]
    y = value[:,1]
    z = value[:,2]
    var = value[:,3]/(0.5*Uinf**2)
    p=pd.DataFrame()
    p["x"]=x
    p["y"]=y
    p["z"]=z
    p["Cp"]=var
    p["p"]=(var*0.5*densityAir*Uinf*Uinf)
    p.to_csv(saveFile,index=False)
    return x, y, z, var, (var*0.5*densityAir*Uinf*Uinf)
    

def runFluid(caseDir, iteration, inBoundaryFiles, pressureFileName, aoa, parallel, 
             Uinf = 51.4815, outBox = (-50,75,-50,50,0.01),endTime=10000, willDelete = False,
             densityAir =1.1839):
    """"
    Args:
        iteration (int) : iteration index
        inBoundaryFiles (list of strings) : boundary nodes provided by a list of in csv file names
        pressureFileName: distribution of pressure over the boundary is saved by this name
        aoa (float): angle of attack
        parallel (boolean): True if the simulation is run in parallel
        Uinf (float): farfield velocity magnitude
        outBox (list, tuple, np.array): to define the farfield boundary, [min-x, max-x, min-y, max-y, thickness following z]
        endTime (int): maximum number of iterations in the simulation
        willDate(boolean): True if the running folder is deleted after finishing
        densityAir (float): density of air to compute the force coefficients
    Returns:
        CL (float), CD (float), maxInt(float): lift coefficient, drag coefficient, the converged iteration 
    """
    boundFiles = [os.path.join(caseDir, ff) for ff in inBoundaryFiles]
    profile = utils.createProfile(boundFiles) 
    if parallel:
        executable="mpiexec -np 8 simpleFoam -parallel"
    else:
        executable = "simpleFoam"
    CL, CD, maxInt = run2DCase.runCase(profile, f"{caseDir}/referenceCase",f"{caseDir}/fluid-run-{iteration}", 
                                   aoa=aoa, outBox = outBox, 
                                   Uinf=Uinf, executable=executable, parallel=parallel,
                                   endTime=endTime, willDelete = willDelete)
    print(f"CL = {CL}, CD = {CD} maxInt = {maxInt}")
    pressureFileName = os.path.join(f"{caseDir}/fluid-run-{iteration}","Cp.csv") 
    x,y, z, Cp, p = getCp(f"{caseDir}/fluid-run-{iteration}/postProcessing/samplePwall1/surface",
                     maxInt,"p_patch_aerofoil.raw",Uinf, densityAir, pressureFileName)
    if not(os.path.exists(pressureFileName)):
        print(f"File {pressureFileName} does not exist")
        return np.nan,np.nan,0
    p = pd.DataFrame()
    p["CL"] = np.array([CL])
    p["CD"] = np.array([CD])
    p["maxInt"] = np.array([maxInt])
    p.to_csv(os.path.join(f"{caseDir}/fluid-run-{iteration}","CLCD.csv"), index=False)
    return CL, CD, maxInt

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    plt.close("all")
    inBoundaryFiles = ["initialGeo/boundaryNodes0.csv",
                 "initialGeo/boundaryNodes1.csv",
                 "initialGeo/boundaryNodes2.csv"]
    CL, CD, maxInt = runFluid(".",1000,inBoundaryFiles,"Cp.csv",10,True) 
    print(f"result = {CL, CD, maxInt}")
                  
    plt.figure()
    profile = utils.createProfile(inBoundaryFiles) 
    plt.plot(profile[:,0],profile[:,1],"g.-", label="deformed airfoil")
    
    aerofoil = Aerofoil("naca0012",*MPXX_family("0012",normalised = True))
    profile2 = aerofoil.getProfile()
    plt.plot(profile2[:,0],profile2[:,1],"r-",label="initial")
    plt.legend()
    
    plt.show()
