from runSolid import *
from runFluid import *
import json

def runFSI(params):
    iteration, maxDisp,caseDir,runFromStart = 0, 0, params["caseDirectory"],params["runFromStart"]
    if runFromStart==False:
        with open("timeLogFile.txt", "a") as f: f.write("\n Starting mid-simulation. Checking entry point.")
        # Find where the simulation has got to by checking the existance of files.
        if os.path.exists(os.path.join(caseDir, "logIter.txt")):
            print("logIter found")
            initialBoundaryFiles = ["initialGeo/"+f for f in os.listdir("initialGeo") if f.startswith("boundary")]
            with open(f"{caseDir}/logIter.txt", "r") as f: # if run hasn't got past first iteration 
                prevData =  f.readlines()[-1].split(",")
                iteration, maxDisp, lastError = int(prevData[0]), float(prevData[1]), float(prevData[2])
            if lastError < params["relativeTolerance"]: return maxDisp
            print(iteration, maxDisp)
            if iteration == -1:
                iteration = 0
                if not os.path.exists(os.path.join(caseDir, f"fluid-run-0", "stress.csv")): 
                    runFromStart=True # first fluid sim hasm't completed
                    with open("timeLogFile.txt", "a") as f: f.write("\n Fluid-run-0 did not complete. Run from start.")
            if (os.path.exists(os.path.join(caseDir, f"fluid-run-{iteration}", "stress.csv")) and 
                not os.path.exists(os.path.join(caseDir, f"solid-run-{iteration}", "extractedData.csv"))):# Fluid part has run but solid has not.
                with open(f"{caseDir}/logIter.txt","a") as logFile:
                    maxDisp, minDistShunt = runSolid(caseDir,iteration,initialBoundaryFiles,
                                        initialBoundaryFiles,params)
                    logFile.write(f"0,{maxDisp},1,{minDistShunt}\n")
                    logFile.flush()
            else: # run the solid part of the first sim
                print(os.system(f"ls {caseDir}"))
                    
        else:
            runFromStart = True
            with open("timeLogFile.txt", "a") as f: f.write("\n No log file present. Run from start.")
    if runFromStart:
        with open("timeLogFile.txt", "a") as f: f.write("\n Running from start. Tidy and remesh initial boundary.")
        os.system(f"rm -rf *.txt {caseDir}; mkdir {caseDir}; python3 makeInitialBoundary.py; \
                   cp -r initialGeo {caseDir}; cp -rf referenceFluid {caseDir};  cp plotGeo.py {caseDir}")
        os.system(f"cp -r referenceSolid {caseDir}")

        initialBoundaryFiles = ["initialGeo/"+f for f in os.listdir("initialGeo") if f.startswith("boundary")]
        with open(f"{caseDir}/logIter.txt","a") as logFile:
            logFile.write(f"Iteration,maxDisp,errorDisp,minDistShunt\n")
            logFile.write(f"-1,0,1,10\n")
            logFile.flush()
            runFluid(caseDir, iteration,initialBoundaryFiles,params)
            maxDisp, minDistShunt = runSolid(caseDir,iteration,initialBoundaryFiles,
                                    initialBoundaryFiles,params)
            logFile.write(f"0,{maxDisp},1,{minDistShunt}\n")
            logFile.flush()
    maxDisp0 = maxDisp

    while True:
        with open("timeLogFile.txt", "a") as f: f.write(f"\n Starting iteration {iteration+1}")
        with open(f"{caseDir}/logIter.txt","a") as logFile: 
            iteration+=1
            inBoundaryFiles = [os.path.join(f"solid-run-{iteration-1}", i) for i in os.listdir(os.path.join(caseDir, f"solid-run-{iteration-1}")) if 'boundary' in i]
            if not os.path.exists(os.path.join(caseDir, f"fluid-run-{iteration}", "stress.csv")):# have the latest fluid but not solid
                runFluid(caseDir, iteration,inBoundaryFiles,params) 
            maxDispPrev = maxDisp
            maxDisp, minDistShunt = runSolid(caseDir,iteration,initialBoundaryFiles,inBoundaryFiles, params)
            errorDisp = np.abs(maxDisp - maxDispPrev)/maxDisp0

            logFile.write(f"{iteration},{maxDisp},{errorDisp},{minDistShunt}\n")
            logFile.flush()
            if errorDisp < params["relativeTolerance"]: 
                print("CONVERGENCE BY ERROR TOL: " + str(params["relativeTolerance"]))
                break
            elif np.abs(maxDisp) < params["absoluteTolerance"]:
                print("CONVERGENCE BY ABSOLUTE TOL: " + str(params["absoluteTolerance"]))
                break
            elif iteration == params["maxIterations"]:
                print("MAX ITERATIONS REACHED")
                break
            elif minDistShunt==0:
                print("COLLISION WITH SHUNT")
                break
        logFile.close()   
    os.system(f"python3 plotGeo.py {caseDir}")
    if params["cleanFiles"]: 
        with open("timeLogFile.txt", "a") as f: f.write(f"\n Removing intermediary files with cleanFiles.py")
        os.system("python3 cleanFiles.py")
    return maxDisp

if __name__ == "__main__":
    # Read in all params from config file
    with open("timeLogFile.txt", "a") as f: f.write("\n Starting FSI simulation :  "+str(datetime.now().strftime("%Y-%m-%d %H-%M-%S")))
    maxDisp = runFSI( json.loads(open("config.json", "r").read()))
    print(f"Result : maxDisp ={maxDisp}")

