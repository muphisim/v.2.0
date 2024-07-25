import os
import json

params = json.loads(open("config.json", "r").read())

def getListMaxInter(folder):
    allDirs = [d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder,d))]
    allNums = []
    key = "solid-run-"
    for d in allDirs:
        if key in d:
            allNums.append(int(d[len(key):]))
    return max(allNums)


for runName in [d for d in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(),d))]:
    runName = os.getcwd()
    folder = os.path.join(runName, params["caseDirectory"])
    print(folder)
    try:
        lastIt = getListMaxInter(folder)
        os.system(f"mv {folder}/fluid-run-{lastIt} {runName}/; rm -rf {folder}/fluid-run-*")
        os.system(f"mv {folder}/solid-run-{lastIt} {runName}/; rm -rf {folder}/solid-run-*")
        os.system(f"cp {folder}/logIter.txt {runName}/;")
    except: 
        print(folder)
        continue
os.system(f"rm -rf {folder}/log-run.txt ")
