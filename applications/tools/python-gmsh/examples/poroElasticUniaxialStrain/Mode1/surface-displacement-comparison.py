import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math

plt.close("all")
plt.figure()

L = 100
Q=33.3
kappa=1e-6
mu = 2.67e-5
nu = 0.35
alpha = 1
G = 7.78
c=1/((mu/kappa)*(((alpha**2)*(1-2*nu))/(2*G*(1-nu))+1/Q))
n = ((alpha*(1-2*nu))/(2*(1-nu)))
a = (1-2*nu)/(2*G*(1-nu))
ai = a/(1+(alpha**2)*a*Q)
K = 21/(3*(1-2*nu))

ps=0.1
pi=ps*(a-ai)/(a)
t = np.linspace(0,10,1000)
y1 = t*0
nlist=np.linspace(1,1000,1000)    

#for i in range(len(t)):
       # for jj in range(len(nlist)):
            #k = np.exp(-c*t[i]*(np.pi*(nlist[jj]*2+1)/(2*L))**2)/((nlist[jj]*2+1)**2) + k 
            #k = (1-np.exp(-nlist[jj]*nlist[jj]*np.pi*np.pi*c*t[i]/(4*L*L)))*8/(nlist[jj]*nlist[jj]*np.pi*np.pi) +k

        #y[i] = -(k*n*ps*L)/G + ai*L*ps 
        #y[i] = 2*(a-ai)*ps*np.sqrt(c*t[i]/np.pi) + ai*L*ps 
        #y[i] = (-8/(np.pi**2))*(a-ai)*L*ps*k + a*L*ps
        #y[i] = L*ps*(1/(K+1.33333333333*G))
        #k = 0

for i in range(len(t)):  
    y1[i] = 2*(a-ai)*ps*np.sqrt(c*t[i]/np.pi) + ai*L*ps #- n*ps*L*4*np.sqrt(c*t[i]/(4*L*L))/G

plt.plot(t,y1,"-",label="Analytical") 


p = pd.read_csv("out/SurfaceDisplacement_Unknown0_Rough.csv") 
x = p["Time"].values 
y = p["Node11786Comp0"].values 
plt.plot(x,y,".",label="MuPhiSim")

y1 = y*0
k = y*0



for i in range(len(y)):
    y1[i] = 2*(a-ai)*ps*np.sqrt(c*t[i]/np.pi) + ai*L*ps #- n*ps*L*4*np.sqrt(c*t[i]/(4*L*L))/G
    if np.abs(y1[i] - y[i]) > 0:
        k[i] = (y1[i] - y[i])/y1[i]

total = sum(np.abs(k))/len(t)
index_max = np.argmax(np.abs(k))
print(max(np.abs(k)))
print(total)
print(index_max) 
print(t[index_max])
print(min(abs(k)))


plt.xlabel("Time, $t$ (s)")
plt.ylabel("Surface displacement, $u_s$ (cm)")  

plt.legend(ncol=1)  

plt.savefig("surface-displacement-coupled.png",dpi=1000)

plt.show()
