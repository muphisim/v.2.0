import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

plt.close("all")
plt.figure()


Q=3.33E-1
kappa=3.75E-5
mu = 1E-5
c=kappa*Q/mu

pi=0
ps=10e6
x = np.linspace(0,100,1000)
y = x*0
tlist=[0.1,1,3,10]
for t in tlist:
    
    print(f"semi-finite creterion {50/np.sqrt(c*t) -2}")
    
    for i in range(len(x)):
      y[i] = ps +(pi-ps)* math.erf(x[i]/2/np.sqrt(c*t))
    
    plt.plot(x,y,"-",label=f"analytic-{t}s")

    
p = pd.read_csv("out-uncoupled/time-0.1.csv") 
x = p["Points_0"].values
y = p["ExtraDof1"].values
plt.plot(x,y,".",label="current-0.1s")

p = pd.read_csv("out-uncoupled/time-1.csv") 
x = p["Points_0"].values
y = p["ExtraDof1"].values
plt.plot(x,y,".",label="current-1s")

p = pd.read_csv("out-uncoupled/time-3.csv") 
x = p["Points_0"].values
y = p["ExtraDof1"].values
plt.plot(x,y,".",label="current-3s")

p = pd.read_csv("out-uncoupled/time-10.csv") 
x = p["Points_0"].values
y = p["ExtraDof1"].values
plt.plot(x,y,".",label="current-10s")

plt.xlabel("$x$-coordinate")
plt.ylabel("p")

plt.legend(ncol=2)

plt.savefig("comparison-with-analytic.png",dpi=1000)

plt.show()
