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
L = 100
c=1/((mu/kappa)*(((alpha**2)*(1-2*nu))/(2*G*(1-nu))+1/Q))
n = alpha*(1-2*nu)/(2*(1-nu))
a = (1-2*nu)/(2*G*(1-nu))
ai = a/(1+(alpha**2)*a*Q)

ps=0.1
pi=ps*(a-ai)/(a)
x = np.linspace(0,100,10000)
y = x*0
u = x*0
tlist=[0.1,1,3,10]
nlist=np.arange(1,1003,2)
print(n)
k = 0
for t in tlist:
    
    print(f"semi-finite criterion {50/np.sqrt(c*t) -2}")
    
    for i in range(len(x)):
    
        for jj in range(len(nlist)):
            k = np.exp(-c*t*(np.pi*nlist[jj]/(2*L))**2)*np.sin(np.pi*nlist[jj]*x[i]/(2*L))/nlist[jj] + k 
            
           
        y[i] = k*(4*ps*(a-ai))/(np.pi*alpha*a)
        #y[i] =  k*(4*(ps+pi)*(a-ai))/(np.pi*alpha*a) 
        #y[i] =  k*(4*(ps-pi))/(np.pi)       
        k = 0
        

    plt.plot(x,y,"-",label=f"analytic-{t}s")


p = pd.read_csv("out/time-0.1s.csv") 
x = p["Points_0"].values 
y = p["ExtraDof1"].values
y1 = y*0
z = y1*0
plt.plot(x,y,".",label="MuPhiSim-0.1s")
for i in range(len(x)):
    y1[i]= ps +(pi-ps)* math.erf(x[i]/2/np.sqrt(c*0.1))
    if y1[i] - y[i] > 0 and y1[i] > 0.00001:
     z[i] = 100*(y1[i] - y[i])/y1[i]
     

total = sum(z)/len(x)
perc = max(z)
index_max = np.argmax(z)
print(total)
print(perc) 
print(index_max) 
print(x[index_max])    

p = pd.read_csv("out/time-1s.csv") 
x = p["Points_0"].values 
y = p["ExtraDof1"].values
plt.plot(x,y,".",label="MuPhiSim-1s")
for i in range(len(x)):
    y1[i]= ps +(pi-ps)* math.erf(x[i]/2/np.sqrt(c*1))
    if y1[i] - y[i] > 0 and y1[i] > 0.00001:
     z[i] = 100*(y1[i] - y[i])/y1[i]
     

total = sum(z)/len(x)
perc = max(z)
index_max = np.argmax(z)
print(total)
print(perc) 
print(index_max) 
print(x[index_max])


p = pd.read_csv("out/time-3s.csv") 
x = p["Points_0"].values 
y = p["ExtraDof1"].values
plt.plot(x,y,".",label="MuPhiSim-3s")
for i in range(len(x)):
    y1[i]= ps +(pi-ps)* math.erf(x[i]/2/np.sqrt(c*3))
    if y1[i] - y[i] > 0 and y1[i] > 0.00001:
     z[i] = 100*(y1[i] - y[i])/y1[i]
     

total = sum(z)/len(x)
perc = max(z)
index_max = np.argmax(z)
print(total)
print(perc) 
print(index_max) 
print(x[index_max])

     


p = pd.read_csv("out/time-10s.csv") 
x = p["Points_0"].values 
y = p["ExtraDof1"].values
plt.plot(x,y,".",label="MuPhiSim-10s")
for i in range(len(x)):
    y1[i]= ps +(pi-ps)* math.erf(x[i]/2/np.sqrt(c*10))
    if y1[i] - y[i] > 0 and y1[i] > 0.00001:
     z[i] = 100*(y1[i] - y[i])/y1[i]
     

total = sum(z)/len(x)
perc = max(z)
index_max = np.argmax(z)
print(total)
print(perc) 
print(index_max) 
print(x[index_max])
    

plt.legend(ncol=2)

#plt.savefig("comparison-coupled-Mode1.png",dpi=1000)

plt.show()
