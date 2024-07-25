import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

plt.close("all")
plt.figure()

p = pd.read_csv("solid-MuPhiSim/output/time-5s.csv") 
x = p["Points_0"].values 
y1 = p["ExtraDof1"].values 
plt.plot(x,y1,".",label="$t$=5s")

p = pd.read_csv("solid-MuPhiSim/output/time-10s.csv") 
x = p["Points_0"].values 
y1 = p["ExtraDof1"].values 
plt.plot(x,y1,".",label="$t$=10s")

p = pd.read_csv("solid-MuPhiSim/output/time-250s.csv") 
x = p["Points_0"].values 
y1 = p["ExtraDof1"].values 
plt.plot(x,y1,".",label="$t$=25s")

p = pd.read_csv("solid-MuPhiSim/output/time-50s.csv") 
x = p["Points_0"].values 
y1 = p["ExtraDof1"].values 
plt.plot(x,y1,".",label="$t$=50s")

p = pd.read_csv("solid-MuPhiSim/output/time-100s.csv") 
x = p["Points_0"].values 
y1 = p["ExtraDof1"].values 
plt.plot(x,y1,".",label="$t$=100s")


plt.xlabel("$x$-coordinate ")
plt.ylabel("Pore pressure, $p$")  

plt.legend(ncol=2)  

plt.show()
