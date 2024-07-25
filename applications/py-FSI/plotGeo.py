import matplotlib.pyplot as plt
import pandas as pd
import random
import os

plt.close("all")

plt.figure()

def getListMaxInter(src):
    allDirs = [d for d in os.listdir(src) if os.path.isdir(os.path.join(src,d))]
    allNums = []
    key = "fluid-run-"
    for d in allDirs:
        if key in d:
            allNums.append(int(d[len(key):]))
    return max(allNums)

numFluidRun=getListMaxInter(os.getcwd())
color = ["#"+''.join([random.choice('0123456789ABCDEF') for j in range(8)])
             for i in range(numFluidRun+2)]


files = ["initialGeo/boundaryNodes0.csv","initialGeo/boundaryNodes1.csv","initialGeo/boundaryNodes2.csv"]

for ff in files:
    p = pd.read_csv(ff)
    label = "initial" if files.index(ff)==0 else ""
    plt.plot(p.values[:,1],p.values[:,2],"-", label=label , color=color[0])

for i in range(numFluidRun+1):
    files = [f"solid-run-{i}/boundaryNodes0.csv",
             f"solid-run-{i}/boundaryNodes1.csv",
             f"solid-run-{i}/boundaryNodes2.csv"]
    for ff in files:
        p = pd.read_csv(ff)
        label = f"solid-run-{i}" if files.index(ff)==0 else ""
        plt.plot(p.values[:,1],p.values[:,2],"-", label=label, color=color[i+1])


plt.legend()
plt.axis("equal")


plt.figure()

for i in range(numFluidRun+1):
    CpFile = f"fluid-run-{i}/Cp.csv"
    Cp = pd.read_csv(CpFile)
    plt.plot(Cp.values[:,0],Cp.values[:,4],".",label=CpFile)



plt.legend()

plt.show()
