import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

plt.close("all")
plt.figure()

p = pd.read_csv("fluid-openFOAM/precice-Fluid-watchpoint-Wall-0.log",delim_whitespace=True) 
x = p["Time"].values 
y1 = p["Pressure0"].values 
plt.plot(x,y1,"-",label="Fluid, Point 1")

p = pd.read_csv("solid-MuPhiSim/precice-Solid-watchpoint-Wall-0.log",delim_whitespace=True) 
x = p["Time"].values 
y2 = p["Pressure0"].values 
plt.plot(x,y2,"--",label="Solid, Point 1")
#y = y2 - y1
#plt.plot(x,y,"-",label="Difference, Point 1")

p = pd.read_csv("fluid-openFOAM/precice-Fluid-watchpoint-Wall-1.log",delim_whitespace=True) 
x = p["Time"].values 
y1 = p["Pressure0"].values 
plt.plot(x,y1,"-",label="Fluid, Point 2")

p = pd.read_csv("solid-MuPhiSim/precice-Solid-watchpoint-Wall-1.log",delim_whitespace=True) 
x = p["Time"].values 
y2 = p["Pressure0"].values 
#y = y2 - y1
plt.plot(x,y2,"--",label="Solid, Point 2")
#plt.plot(x,y,"-",label="Difference, Point 2")

p = pd.read_csv("fluid-openFOAM/precice-Fluid-watchpoint-Wall-2.log",delim_whitespace=True) 
x = p["Time"].values 
y1 = p["Pressure0"].values 
plt.plot(x,y1,"-",label="Fluid, Point 3")

p = pd.read_csv("solid-MuPhiSim/precice-Solid-watchpoint-Wall-2.log",delim_whitespace=True) 
x = p["Time"].values 
y2 = p["Pressure0"].values 
#y = y2 - y1
plt.plot(x,y2,"--",label="Solid, Point 3")
#plt.plot(x,y,"-",label="Difference, Point 3")

p = pd.read_csv("fluid-openFOAM/precice-Fluid-watchpoint-Wall-3.log",delim_whitespace=True) 
x = p["Time"].values 
y1 = p["Pressure0"].values 
plt.plot(x,y1,"-",label="Fluid, Point 4")

p = pd.read_csv("solid-MuPhiSim/precice-Solid-watchpoint-Wall-3.log",delim_whitespace=True) 
x = p["Time"].values 
y2 = p["Pressure0"].values 
#y = y2 - y1
plt.plot(x,y2,"--",label="Solid, Point 4")
#plt.plot(x,y,"-",label="Difference, Point 4")

plt.xlabel("Time, $t$ (s)")
plt.ylabel("Pressure in $x$")  

plt.legend(ncol=2)  

plt.show()
