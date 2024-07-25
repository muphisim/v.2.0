import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

plt.close("all")
plt.figure()


Q=33.3
kappa=1e-6
mu = 2.67e-5
nu = 0.35
alpha = 1
G = 7.78

c=1/((mu/kappa)*(((alpha**2)*(1-2*nu))/(2*G*(1-nu))+1/Q))

pi=0
ps=0.1
x = np.linspace(0,100,100000)
y = x*0
tlist=[0.1,1,3,10]
for t in tlist:
    
    print(f"semi-finite criterion {50/np.sqrt(c*t) -2}")
    
    for i in range(len(x)):
      y[i] = ps +(pi-ps)* math.erf(x[i]/2/np.sqrt(c*t))
    
    plt.plot(x,y,"-",label=f"analytic-{t}s")

    
p = pd.read_csv("out/time-0.1s.csv") 
x = p["Points_0"].values
y = p["ExtraDof1"].values
plt.plot(x,y,".",label="MuPhiSim-0.1s")

p = pd.read_csv("out/time-1s.csv") 
x = p["Points_0"].values
y = p["ExtraDof1"].values
plt.plot(x,y,".",label="MuPhiSim-1s")

p = pd.read_csv("out/time-3s.csv") 
x = p["Points_0"].values
y = p["ExtraDof1"].values
plt.plot(x,y,".",label="MuPhiSim-3s")

p = pd.read_csv("out/time-10s.csv") 
x = p["Points_0"].values
y = p["ExtraDof1"].values
plt.plot(x,y,".",label="MuPhiSim-10s")

plt.xlabel("$x$-coordinate")
plt.ylabel("Pore pressure variation, $p$")

plt.legend(ncol=2)

plt.savefig("comparison-coupled-Mode2.png",dpi=1000)

plt.show()
