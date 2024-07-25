import numpy as np

def Area(profile):
    """ compute area using Green's theorem
    """
    A = 0
    pFull = np.append(profile,[profile[0].tolist()], axis=0)
    numSegs = len(pFull)-1
    for i in range(numSegs):
        xc = (pFull[i,0]+pFull[i+1,0])/2.
        dx = pFull[i+1,0]-pFull[i,0]
        yc = (pFull[i,1]+pFull[i+1,1])/2.
        dy = pFull[i+1,1]-pFull[i,1]
        A += 0.5*(xc*dy - yc*dx)
    return A
    
if __name__ == "__main__":
    import matplotlib.pyplot as plt    
    
    theta =np.linspace(0,2*np.pi,200, endpoint=False)
    rx = 1
    ry = 0.2
    x = rx*np.cos(theta)
    y = ry*np.sin(theta)
    profile = np.stack((x,y),axis=1)
    a = Area(profile)
    print(f"approx={a} true value = {np.pi*rx*ry}")
    
    profile =  np.stack((x,y),axis=1)

    plt.plot(profile[:,0],profile[:,1],"r--")
    plt.axis("equal")
    plt.show()
