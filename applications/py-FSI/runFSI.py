from runSolid import *
from runFluid import *

def runFSI(caseDir, initialBoundaryFiles, aoa, parallel, E, nu, relaxFactor, maxIter, tol, absTol):
    print(locals())
    logFile = open(f"{caseDir}/logIter.txt","w")
    logFile.write(f"Iteration,maxDisp,errorDisp,maxForce,errorForce,CL,errorCL,CD,errorCD\n")
    logFile.flush()
    iteration=0
    CL, CD, maxInt = runFluid(caseDir, iteration,initialBoundaryFiles,"Cp.csv",aoa,parallel) 
    CpFile = f"{caseDir}/fluid-run-0/Cp.csv"
    gmshModel, pressureData, maxDisp, maxForce = runSolid(caseDir,iteration,initialBoundaryFiles,CpFile,[(1,4),(1,5)],[E,nu],relaxFactor)
    logFile.write(f"0,{maxDisp},1,{maxForce},1,{CL},1,{CD},1\n")
    logFile.flush()
    
    if np.abs(maxDisp) < absTol:
        print("converge by absolute tol of maxDisp")
        logFile.close()  
        return CL, CD, maxInt 

    
    maxDisp0 = maxDisp
    maxForce0 = maxForce
    CL0=CL
    CD0=CD
    while True:
        iteration+=1
        inBoundaryFiles = [f"solid-run-{iteration-1}/boundaryNodes0.csv",
                     f"solid-run-{iteration-1}/boundaryNodes1.csv",
                     f"solid-run-{iteration-1}/boundaryNodes2.csv"]
        CLPrev=CL
        CDPrev=CD
        CL, CD, maxInt = runFluid(caseDir, iteration,inBoundaryFiles,"Cp.csv",aoa,parallel) 
        CpFile = f"{caseDir}/fluid-run-{iteration}/Cp.csv"
        maxDispPrev = maxDisp
        maxForcePrev = maxForce
        gmshModel, pressureData, maxDisp, maxForce = runSolid(caseDir,iteration,initialBoundaryFiles,CpFile,[(1,4),(1,5)],[E,nu],relaxFactor)
        #
        errorDisp = np.abs(maxDisp - maxDispPrev)/maxDisp0
        errorForce = np.abs(maxForce - maxForcePrev)/maxForce0
        errorCL=np.abs(CL-CLPrev)/CL0
        errorCD=np.abs(CD-CDPrev)/CD0
        logFile.write(f"{iteration},{maxDisp},{errorDisp},{maxForce},{errorForce},{CL},{errorCL},{CD},{errorCD}\n")
        logFile.flush()
        if errorDisp < tol: 
            print("convergence by errorDisp!!!")
            break
        elif (errorCL < tol) and (errorCD < tol):
            print("convergence by errorCL and errorCD!!!")
            break
        elif np.abs(maxDisp) < absTol:
            print("converge by absolute tol of maxDisp!!!")
            break
        elif iteration == maxIter:
            print("maximal iterations is reached!!!")
            break
    logFile.close()   
    return CL, CD, maxInt

if __name__ == "__main__":

    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-case", "--case", type=str, required=True,
                        help="case dir", metavar="caseDir")
    parser.add_argument("-aoa", "--aoa", type=float, required=True,
                        help="angle of attach", metavar="aoa")
    parser.add_argument("-parallel", "--parallel", type=int, required=True,
                        help="parallel", metavar="0/1") 
    parser.add_argument("-young", "--young", type=float, required=True,
                        help="young modulus", metavar="E")  
    parser.add_argument("-nu", "--nu", type=float, required=True,
                        help="Poisson ratio", metavar="nu")        
    parser.add_argument("-omega", "--omega", type=float, required=True,
                        help="relax factor", metavar="omega")  
    parser.add_argument("-maxIter", "--maxIter", type=int, required=False,
                        help="max iteration", metavar="maxIter")  
    parser.add_argument("-tol", "--tol", type=int, required=False,
                        help="relative tolerance", metavar="tol")      
    parser.add_argument("-absTol", "--absTol", type=int, required=False,
                        help="absolute tolerance", metavar="absTol")
                                             
    args = parser.parse_args()
    print(args)
    aoa=args.aoa
    parallel= bool(args.parallel)
    E = args.young
    nu = args.nu
    omega = args.omega
    caseDir = args.case
    
    maxIter = 10 
    tol = 5e-3
    absTol = 5e-4
    
    if args.maxIter is not None:
        maxIter = args.maxIter
    if args.tol is not None:
        tol=args.tol
    if args.absTol is not None:
        absTol = args.absTol
    
    
    initialBoundaryFiles = [f"initialGeo/boundaryNodes0.csv",
                 "initialGeo/boundaryNodes1.csv",
                 "initialGeo/boundaryNodes2.csv"]  
               
    CL, CD, maxInt = runFSI(caseDir, initialBoundaryFiles,aoa,parallel,E,nu, omega,maxIter, tol, absTol)
    print(f"Result : CL ={CL}, CD={CD}, maxInt ={maxInt}")
    p = pd.DataFrame()
    p["CL"]=[CL]
    p["CD"]=[CD]
    p.to_csv(f"{caseDir}/finalCLCD.csv",index=False)

