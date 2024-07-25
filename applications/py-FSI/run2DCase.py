import numpy as np
import os
import time
import subprocess

def writeGeoFile_fine(profile, geoFileName, aoa, xmin, xmax, zmin, zmax, thickness, algo):
    """write geometry file with gmsh 
        Args:
            geoFileName (string)    : name of geometry file
            xmin, xmax, zmin, zmax, thickness (floats): define the out-box for farfield condition
            
            algo: msh algorithm, 1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 11: Quasi-structured Quad
            
        Returns:
            None
            
    """
    x= profile[:,0]
    z = profile[:,1]
    npoints = len(x)
    file = open(geoFileName,"w")
    file.write("General.ExpertMode = 1;\n")
    file.write(f"Mesh.Algorithm = {algo};\nMesh.RecombinationAlgorithm = 0;\nMesh.RecombineAll = 1;\n")
    file.write("lsca=0.0015;\n")
    for ii in range(npoints):
        file.write(f"Point({ii+1})={{{x[ii]},{z[ii]},{-0.5*thickness},lsca}};\n")
    
    posUp = -1
    posDown = -1
    for ii in range(npoints):
        if x[ii] < 0.5 and posUp<0:
            posUp = ii
        if posUp >0  and x[ii] >0.5:
            posDown = ii
        if posUp >0 and posDown>0:
            break
    print(f"posUp={posUp}, posDown={posDown} npoints = {npoints}")
        
    file.write(f"ls1 = newl; Spline(ls1)={{1:{posUp}}};\n")
    file.write(f"ls2 = newl; Spline(ls2)={{{posUp}:{posDown}}};\n")
    file.write(f"ls3 = newl; Spline(ls3)={{{posDown}:{npoints}, 1}};\n")
    file.write("Transfinite Curve {ls1, -ls3} = 101 Using Progression 1.025;\n")
    file.write("Transfinite Curve {ls2} = 301 Using Bump 20.04;\n")
    file.write("lsca2 = 1.5;\n")
    file.write(f"p1 = newp; Point(p1) = {{{xmin},  {zmin}, {-0.5*thickness},lsca2}};\n")
    file.write(f"p2 = newp; Point(p2) = {{{xmax},  {zmin}, {-0.5*thickness},lsca2}};\n")
    file.write(f"p3 = newp; Point(p3) = {{{xmax},  {zmax}, {-0.5*thickness},lsca2}};\n")
    file.write(f"p4 = newp; Point(p4) = {{{xmin},  {zmax}, {-0.5*thickness},lsca2}};\n")
    file.write(f"p5 = newp; Point(p5) = {{{xmax},  {xmax*np.sin(aoa*np.pi/180.)},{-0.5*thickness}, 0.1*lsca2}};\n")
    file.write(f"p6 = newp; Point(p6) = {{{xmin},  {xmin*np.sin(aoa*np.pi/180.)},{-0.5*thickness}, 0.1*lsca2}};\n")
    file.write("l1=newl; Line(l1)={p1,p2};\n")
    file.write("l2=newl; Line(l2)={p2,p5};\n")
    file.write("l3=newl; Line(l3)={p5,p3};\n")
    file.write("l4=newl; Line(l4)={p3,p4};\n")
    file.write("l5=newl; Line(l5)={p4,p6};\n")
    file.write("l6=newl; Line(l6)={p6,p1};\n")
    file.write("loop1=newl;\nLine Loop(loop1)={ls1,ls2, ls3};\n")
    file.write("loop2=newl;\nLine Loop(loop2)={l1,l2,l3,l4,l5,l6};\n")
    file.write("s1 = news; Plane Surface(s1)={loop1,loop2};\n")
    file.write(f"out[] = Extrude {{0, 0, {thickness}}} {{ Surface{{s1}}; Layers {{1}}; Recombine;}};\n")
    file.write("Physical Volume(1) = {out[1]};\n")
    file.write(f"Physical Surface(\"front\",1) = {{s1}};\n")
    file.write("Physical Surface(\"back\",2) = {out[0]};\n")
    file.write("Physical Surface(\"aerofoil\",3) = {out[2],out[3],out[4]};\n")
    file.write("Physical Surface(\"inlet\",4) = {out[5],out[6],out[7],out[10]};\n")
    file.write("Physical Surface(\"outlet\",5) = {out[8],out[9]};\n")
    file.write("""//Define Boundary Layer
Field[1] = BoundaryLayer;
Field[1].EdgesList = {ls1,ls2,ls3};
Field[1].AnisoMax = 1.0;
Field[1].FanNodesList = {1};
Field[1].hfar = 0.001;
Field[1].hwall_n = 0.0001;
Field[1].thickness = 0.1;
Field[1].Quads = 1;
Field[1].IntersectMetrics = 0;
BoundaryLayer Field = 1;
Mesh.BoundaryLayerFanPoints = 50;""")
    file.close()
    
def writeGeoFile_Bis(profile, geoFileName, aoa, xmin, xmax, zmin, zmax, thickness, algo):
    """write geometry file with gmsh 
        Args:
            geoFileName (string)    : name of geometry file
            xmin, xmax, zmin, zmax, thickness (floats): define the out-box for farfield condition
            
            algo: msh algorithm, 1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 11: Quasi-structured Quad
            
        Returns:
            None
            
    """
    x= profile[:,0]
    z = profile[:,1]
    npoints = len(x)
    file = open(geoFileName,"w")
    file.write("General.ExpertMode = 1;\n")
    file.write(f"Mesh.Algorithm = {algo};\nMesh.RecombinationAlgorithm = 0;\nMesh.RecombineAll = 1;\n")
    file.write("lsca=0.0015;\n")
    for ii in range(npoints):
        file.write(f"Point({ii+1})={{{x[ii]},{z[ii]},{-0.5*thickness},lsca}};\n")
    
    posUp = -1
    posDown = -1
    for ii in range(npoints):
        if x[ii] < 0.5 and posUp<0:
            posUp = ii
        if posUp >0  and x[ii] >0.5:
            posDown = ii
        if posUp >0 and posDown>0:
            break
    print(f"posUp={posUp}, posDown={posDown} npoints = {npoints}")
        
    file.write(f"ls1 = newl; Spline(ls1)={{1:{posUp}}};\n")
    file.write(f"ls2 = newl; Spline(ls2)={{{posUp}:{posDown}}};\n")
    file.write(f"ls3 = newl; Spline(ls3)={{{posDown}:{npoints}, 1}};\n")
    file.write("Transfinite Curve {ls1, -ls3} = 101 Using Progression 1.025;\n")
    file.write("Transfinite Curve {ls2} = 301 Using Bump 20.04;\n")
    file.write("lsca2 = 1.5;\n")
    file.write(f"p1 = newp; Point(p1) = {{{xmin},  {zmin}, {-0.5*thickness},lsca2}};\n")
    file.write(f"p2 = newp; Point(p2) = {{{xmax},  {zmin}, {-0.5*thickness},lsca2}};\n")
    file.write(f"p3 = newp; Point(p3) = {{{xmax},  {zmax}, {-0.5*thickness},lsca2}};\n")
    file.write(f"p4 = newp; Point(p4) = {{{xmin},  {zmax}, {-0.5*thickness},lsca2}};\n")
    file.write(f"p5 = newp; Point(p5) = {{{xmax},  {xmax*np.sin(aoa*np.pi/180.)},{-0.5*thickness}, 0.1*lsca2}};\n")
    file.write(f"p6 = newp; Point(p6) = {{{xmin},  {xmin*np.sin(aoa*np.pi/180.)},{-0.5*thickness}, 0.1*lsca2}};\n")
    file.write("l1=newl; Line(l1)={p1,p2};\n")
    file.write("l2=newl; Line(l2)={p2,p5};\n")
    file.write("l3=newl; Line(l3)={p5,p3};\n")
    file.write("l4=newl; Line(l4)={p3,p4};\n")
    file.write("l5=newl; Line(l5)={p4,p6};\n")
    file.write("l6=newl; Line(l6)={p6,p1};\n")
    file.write("loop1=newl;\nLine Loop(loop1)={ls1,ls2, ls3};\n")
    file.write("loop2=newl;\nLine Loop(loop2)={l1,l2,l3,l4,l5,l6};\n")
    file.write("s1 = news; Plane Surface(s1)={loop1,loop2};\n")
    file.write(f"out[] = Extrude {{0, 0, {thickness}}} {{ Surface{{s1}}; Layers {{1}}; Recombine;}};\n")
    file.write("Physical Volume(1) = {out[1]};\n")
    file.write(f"Physical Surface(\"front\",1) = {{s1}};\n")
    file.write("Physical Surface(\"back\",2) = {out[0]};\n")
    file.write("Physical Surface(\"aerofoil\",3) = {out[2],out[3],out[4]};\n")
    file.write("Physical Surface(\"inlet\",4) = {out[5],out[6],out[7],out[10]};\n")
    file.write("Physical Surface(\"outlet\",5) = {out[8],out[9]};\n")
    file.write("""//Define Boundary Layer
Field[1] = BoundaryLayer;
Field[1].EdgesList = {ls1,ls2, ls3};
Field[1].AnisoMax = 1.0;
Field[1].FanNodesList = {1};
Field[1].hfar = 0.01;
Field[1].hwall_n = 0.0005;
Field[1].thickness = 0.03;
Field[1].Quads = 1;
Field[1].IntersectMetrics = 0;
BoundaryLayer Field = 1;
Mesh.BoundaryLayerFanPoints = 30;""")
    file.close()

def writeGeoFile(profile, geoFileName, aoa, xmin, xmax, zmin, zmax, thickness, algo):
    """write geometry file with gmsh 
        Args:
            geoFileName (string)    : name of geometry file
            xmin, xmax, zmin, zmax, thickness (floats): define the out-box for farfield condition
            
            algo: msh algorithm, 1: MeshAdapt, 2: Automatic, 3: Initial mesh only, 5: Delaunay, 6: Frontal-Delaunay, 7: BAMG, 8: Frontal-Delaunay for Quads, 9: Packing of Parallelograms, 11: Quasi-structured Quad
            
        Returns:
            None
            
    """
    x= profile[:,0]
    z = profile[:,1]
    npoints = len(x)
    file = open(geoFileName,"w")
    file.write("General.ExpertMode = 1;\n")
    file.write(f"Mesh.Algorithm = {algo};\nMesh.RecombinationAlgorithm = 0;\nMesh.RecombineAll = 1;\n")
    file.write("lsca=0.0015;\n")
    for ii in range(npoints):
        file.write(f"Point({ii+1})={{{x[ii]},{z[ii]},{-0.5*thickness},lsca}};\n")
    file.write(f"ls1 = newl; Spline(ls1)={{1:{npoints//2}}};\n ls2 = newl; Spline(ls2)={{{npoints//2}:{npoints}, 1}};\n")
    file.write("Transfinite Curve {ls1, ls2} = 101 Using Bump 0.04;\n")
    file.write("lsca2 = 1.5;\n")
    file.write(f"p1 = newp; Point(p1) = {{{xmin},  {zmin}, {-0.5*thickness},lsca2}};\n")
    file.write(f"p2 = newp; Point(p2) = {{{xmax},  {zmin}, {-0.5*thickness},lsca2}};\n")
    file.write(f"p3 = newp; Point(p3) = {{{xmax},  {zmax}, {-0.5*thickness},lsca2}};\n")
    file.write(f"p4 = newp; Point(p4) = {{{xmin},  {zmax}, {-0.5*thickness},lsca2}};\n")
    file.write(f"p5 = newp; Point(p5) = {{{xmax},  {xmax*np.sin(aoa*np.pi/180.)},{-0.5*thickness}, 0.1*lsca2}};\n")
    file.write(f"p6 = newp; Point(p6) = {{{xmin},  {xmin*np.sin(aoa*np.pi/180.)},{-0.5*thickness}, 0.1*lsca2}};\n")
    file.write("l1=newl; Line(l1)={p1,p2};\n")
    file.write("l2=newl; Line(l2)={p2,p5};\n")
    file.write("l3=newl; Line(l3)={p5,p3};\n")
    file.write("l4=newl; Line(l4)={p3,p4};\n")
    file.write("l5=newl; Line(l5)={p4,p6};\n")
    file.write("l6=newl; Line(l6)={p6,p1};\n")
    file.write("loop1=newl;\nLine Loop(loop1)={ls1,ls2};\n")
    file.write("loop2=newl;\nLine Loop(loop2)={l1,l2,l3,l4,l5,l6};\n")
    file.write("s1 = news; Plane Surface(s1)={loop1,loop2};\n")
    file.write(f"out[] = Extrude {{0, 0, {thickness}}} {{ Surface{{s1}}; Layers {{1}}; Recombine;}};\n")
    file.write("Physical Volume(1) = {out[1]};\n")
    file.write(f"Physical Surface(\"front\",1) = {{s1}};\n")
    file.write("Physical Surface(\"back\",2) = {out[0]};\n")
    file.write("Physical Surface(\"aerofoil\",3) = {out[2],out[3]};\n")
    file.write("Physical Surface(\"inlet\",4) = {out[4],out[5],out[6],out[9]};\n")
    file.write("Physical Surface(\"outlet\",5) = {out[7],out[8]};\n")
    file.write("""//Define Boundary Layer
Field[1] = BoundaryLayer;
Field[1].EdgesList = {ls1,ls2};
Field[1].AnisoMax = 1.0;
Field[1].FanNodesList = {1};
Field[1].hfar = 0.01;
Field[1].hwall_n = 0.0005;
Field[1].thickness = 0.03;
Field[1].Quads = 1;
Field[1].IntersectMetrics = 0;
BoundaryLayer Field = 1;
Mesh.BoundaryLayerFanPoints = 30;""")
    file.close()
    
def runCase(profile, refCase, caseName, outBox, aoa, Uinf, executable, endTime, parallel=False, willDelete = True, fineMesh=False):
    """Run a simulation
        Args:
            refCase : A folder containing a complete reference case, possible without mesh
            caseName: the folder name of the run-case
            outBox:  a list defines far-field range (xmin, xmax, ymin, ymax, thickness)
            aoa     : angle of attach
            Unif    : magnitude of the far-field velocity 
            executable: OpenFOAM executable used to run, eg. simpleFoam
            endTime :end time 
        Returns:
            None
    """

    Ux = Uinf*np.cos(aoa*np.pi/180.)
    Uz = Uinf*np.sin(aoa*np.pi/180.)
    nLift =  [-np.sin(aoa*np.pi/180.), np.cos(aoa*np.pi/180.)]
    nDrag = [np.cos(aoa*np.pi/180.), np.sin(aoa*np.pi/180.)]
    os.system(f"rm -rf {caseName}")
    os.system(f"cp -rf {refCase} {caseName}")
    geoFile = f"{caseName}/aerofoil.geo"
    mshFile = f"{caseName}/aerofoil.msh"
    if fineMesh:
        writeGeoFile_fine(profile,geoFile,aoa, *outBox, algo=8)
    else:   
        writeGeoFile_Bis(profile,geoFile,aoa, *outBox, algo=8)
    
    try:
        outProcess = subprocess.run(["gmsh","-3", f"{geoFile}", f"{mshFile}"],timeout=60) # time out 600p
        ok = outProcess.returncode
    except subprocess.TimeoutExpired:
        print('process ran too long')
        ok = 10000

    if  ok != 0:
        writeGeoFile(profile,geoFile,aoa, *outBox, algo=6)
        try:
            outProcess = subprocess.run(["gmsh","-3", f"{geoFile}", f"{mshFile}"],timeout=60) # time out 600p
            ok = outProcess.returncode
        except subprocess.TimeoutExpired:
            print('process ran too long')
            ok = 10000
        if ok != 0:
            print("cannot mesh")
            os.system(f"rm -rf {caseName}")
            return np.nan,np.nan,0
            
    ok = os.system(f"gmshToFoam {mshFile} -case {caseName}")    
    if ok != 0:
        print("Could not call gmshToFoam")
        os.system(f"rm -rf {caseName}")
        return np.nan,np.nan,0
    os.system(f"checkMesh -case {caseName}")
    #modify U
    try:
        Ufile = open(f"{caseName}/0/U","r")
    except OSError:
        print("Could not open/read file")
        os.system(f"rm -rf {caseName}")
        return np.nan,np.nan,0
    allLines = Ufile.readlines()
    Ufile.close()
    lineIndex = -1
    while True:
        lineIndex += 1
        if lineIndex > len(allLines)-1:
            break
        if ("internalField" in allLines[lineIndex]) :
            allLines[lineIndex] = f"internalField   uniform ({Ux} {Uz} 0.00);\n"
            break
    Ufile = open(f"{caseName}/0/U","w")
    Ufile.writelines(allLines)
    Ufile.close()
    #modify controlDict
    try:
        controlDictfile = open(f"{caseName}/system/controlDict","r")
    except OSError:
        print("Could not open/read file")
        return np.nan,np.nan,0
    allLines = controlDictfile.readlines()
    controlDictfile.close()
    lineIndex = -1
    while True:
        lineIndex += 1
        if lineIndex > len(allLines)-1:
            break
        if "endTime" in allLines[lineIndex]:
            allLines[lineIndex] = f"endTime         {endTime};"
            lineIndex += 1
        if ("forceCoeffs1" in allLines[lineIndex]) :
            allLines[lineIndex+10] = f"        liftDir         ({nLift[0]} {nLift[1]} 0);\n"
            allLines[lineIndex+11] = f"        dragDir         ({nDrag[0]} {nDrag[1]} 0);\n"
            allLines[lineIndex+14] = f"        magUInf         {Uinf};\n"
            allLines[lineIndex+16] = f"        Aref            {outBox[-1]};\n"
            break
    controlDictfile = open(f"{caseName}/system/controlDict","w")
    controlDictfile.writelines(allLines)
    controlDictfile.close()
    # modify bc
    try:
        boundaryFile = open(f"{caseName}/constant/polyMesh/boundary","r")
    except OSError:
        print("Could not open/read file")
        return np.nan,np.nan,0
    allLines = boundaryFile.readlines()
    boundaryFile.close()
    lineIndex = -1
    while True:
        lineIndex += 1
        if lineIndex > len(allLines)-1:
            break
        if ("front" in allLines[lineIndex]) or ("back" in allLines[lineIndex]):
            allLines[lineIndex+2] = allLines[lineIndex+2].replace("patch","empty")
            allLines[lineIndex+3] = allLines[lineIndex+3].replace("patch","empty")
            lineIndex += 3
        
        if "aerofoil" in allLines[lineIndex]:
            allLines[lineIndex+2] = allLines[lineIndex+2].replace("patch","wall")
            allLines[lineIndex+3] = allLines[lineIndex+3].replace("patch","wall")
            lineIndex += 3# -*- coding: utf-8 -*-
    boundaryFile = open(f"{caseName}/constant/polyMesh/boundary","w")
    boundaryFile.writelines(allLines)
    boundaryFile.close()
    print(f"start running {caseName}")
    if parallel:
        os.system(f"decomposePar -case {caseName} >{caseName}/log-decomposePar")
        os.system(f"{executable} -case {caseName} >{caseName}/log-run")
    else:
        os.system(f"{executable} -case {caseName} >{caseName}/log-{executable}")
    print(f"done running {caseName}")
    
    time.sleep(0.1)
    try:
        file = open(f"{caseName}/postProcessing/forceCoeffs1/0/coefficient.dat")
    except OSError:
        print("Could not open/read file")
        os.system(f"rm -rf {caseName}")
        return np.nan,np.nan,0
    lines = file.readlines()
    locMaxInt = None
    locCL = None
    locCD = None 
    for ll in lines:
        if ("# Time" in ll) and ("Cl" in ll) and ("Cd" in ll):
            print(ll)
            res_str = ll.replace('#', '') 
            columns = res_str.split() 
            print(columns)
            locMaxInt = columns.index("Time")
            locCL = columns.index("Cl")
            locCD = columns.index("Cd")
            break
    if (locMaxInt is None) or (locCD is None) or (locCL is None):
        return np.nan,np.nan,0
    
    lastLine = lines[-1].split()
    print(lastLine)
    maxInt= int(lastLine[locMaxInt])
    CD = float(lastLine[locCD])
    CL = float(lastLine[locCL])
    file.close()
    print(f"done: CL = {CL}, CD = {CD} maxInt = {maxInt}")
    if willDelete:
        os.system(f"rm -rf {caseName}")
    if abs(CL)>5 or abs(CD)>5:
        # wrong result
        return np.nan,np.nan,0
    else:
        return CL, CD, maxInt
