import numpy as np
from scipy import interpolate
from scipy import optimize
import Area

class Aerofoil():
    def __init__(self, name, profileUpper, profileLower):
        """
        Args:
            name (string) : the name of profile
            profileUpper :a two dimensional array [[x0,y0], ...] assending order of the first column
            profileLower :a two dimensional array [[x0,y0], ...] assending order of the first column
        """
        self.Type = name
        self.xu = np.array(profileUpper[:,0])
        self.yu = np.array(profileUpper[:,1])
        self.xl = np.array(profileLower[:,0])
        self.yl = np.array(profileLower[:,1])
        self.profile = np.stack((np.append(np.flip(self.xu),self.xl[1:-1]), np.append(np.flip(self.yu),self.yl[1:-1])),axis=1)
        self.interpolaionU = interpolate.splprep([self.xu, self.yu],s=0)
        self.interpolaionL = interpolate.splprep([self.xl, self.yl],s=0)
        
    
    def getArea(self):
        """ compute area using Green's theorem
        """
        return Area.Area(self.profile)
    
    def getProfileUpper(self):
        return np.stack((self.xu, self.yu),axis=1)
    
    def getProfileLower(self):
        return np.stack((self.xl, self.yl),axis=1)
    
        
    def update(self, profileUpper, profileLower):
        """update profile"""
        print("update profile")
        self.xu = np.array(profileUpper[:,0])
        self.yu = np.array(profileUpper[:,1])
        self.xl = np.array(profileLower[:,0])
        self.yl = np.array(profileLower[:,1])
        self.profile = np.stack((np.append(np.flip(self.xu),self.xl[1:-1]), np.append(np.flip(self.yu),self.yl[1:-1])),axis=1)
        self.interpolaionU = interpolate.splprep([self.xu, self.yu],s=0)
        self.interpolaionL = interpolate.splprep([self.xl, self.yl],s=0)
                
    def minDistance(self, p0):
        def val(t, tck):
            x, y = interpolate.splev(t,tck)
            return np.sqrt((x[0]-p0[0])**2 + (y[0]-p0[1])**2)
        
        def func(t, tck):
            x, y = interpolate.splev(t,tck)
            dx, dy = interpolate.splev(t,tck, der=1)
            return ((x[0]-p0[0])*dx[0] +  (y[0]-p0[1])*dy[0])
        
        minUp = np.array(np.sqrt((self.xu-p0[0])**2 + (self.yu-p0[1])**2))
        minLow = np.array(np.sqrt((self.xl-p0[0])**2 + (self.yl-p0[1])**2))
        posMinU, posMinL = np.argmin(minUp), np.argmin(minLow)
        tMinU, tMinL = self.interpolaionU[1][posMinU], self.interpolaionL[1][posMinL]
        valMinUp = minUp[posMinU]
        valMinLow = minLow[posMinL]
        #print("pre:", valMinLow, valMinUp)
        if (tMinU >0) and (tMinU < 1):
            resUp = optimize.root(func, tMinU, args=self.interpolaionU[0])
            valMinUp = val(resUp.x, self.interpolaionU[0])
        
        if (tMinL >0) and (tMinL < 1):
            resLow = optimize.root(func, tMinL, args=self.interpolaionL[0])
            valMinLow = val(resLow.x, self.interpolaionL[0])
        #print("after:", valMinLow, valMinUp)
        return min(valMinLow, valMinUp) 
        
    def getPoint(self, t):
        xu, yu = interpolate.splev(t,self.interpolaionU[0])
        xl, yl =  interpolate.splev(t,self.interpolaionL[0])
        return np.stack((np.append(np.flip(xu),xl), np.append(np.flip(yu),yl)),axis=1)
    
    def diff_list(self, points):
        npoints = len(points)
        res = 0
        for i in range(npoints):
            res +=self.minDistance(points[i])
        return res
    
    def diff(self, other_aerofoil):
        return self.diff_list(other_aerofoil.getProfile())
        
    def getProfile(self):
        """return the profile
        """
        return self.profile
    
    def getProfile_part(self, x_min, x_max, nPoints):
        """return the upper and lower profiles from xc_min to xc_max, nPoints is the sampling points
        """
        tmin = self.__findTFromX(self.interpolaionU[1],self.xu, x_min)
        tmax = self.__findTFromX(self.interpolaionU[1],self.xu, x_max)
        theta = np.linspace(0,np.pi,nPoints)
        t = tmin[0] + 0.5*(1-np.cos(theta))*(tmax[0]-tmin[0])
        xUp, yUp =  interpolate.splev(t,self.interpolaionU[0])

        tmin = self.__findTFromX(self.interpolaionL[1],self.xl, x_min)
        tmax = self.__findTFromX(self.interpolaionL[1],self.xl, x_max)
        t = tmin[0] + 0.5*(1-np.cos(theta))*(tmax[0]-tmin[0])
        xLow, yLow = interpolate.splev(t,self.interpolaionL[0])

        return xUp, yUp, xLow, yLow
    
    @staticmethod
    def __findTFromX(tList, xList, x):
        t = []
        sizeT = len(tList)
        for i in range(sizeT-1):
            if (xList[i] <= x) and (x <= xList[i+1]):
                t.append(tList[i]+ (tList[i+1]-tList[i])*(x-xList[i])/(xList[i+1]-xList[i]))
        return t
            
        

def MPXX_family(MPXX,  xc=None, normalised = True):
    """
        This NACA airfoil series is controlled by 4 digits e.g. NACA 2412, which designate the camber, position of the maximum camber and thickness. If an airfoil number is
        NACA MPXX
        Args:
            M is the maximum camber divided by 100. In the example M=2 so the camber is 0.02 or 2% of the chord
            P is the position of the maximum camber divided by 10. In the example P=4 so the maximum camber is at 0.4 or 40% of the chord.
            XX is the thickness divided by 100. In the example XX=12 so the thiickness is 0.12 or 12% of the chord.
            
            normalised: a flag to make a length of airfoil=1
        Returns:
            xu, yu, xl, yl
        
    """
    m = float(MPXX[0])/100.
    p = float(MPXX[1])/10
    t = float(MPXX[2:])/100
    
    if xc is None:
        theta = np.linspace(0,np.pi,101)
        xc = 0.5*(1-np.cos(theta))
    
    yc = np.array(xc)*0
    dycdxc = np.array(xc)*0
    yt = np.array(xc)*0
    for i in range(len(xc)):
        if xc[i] < p:
            yc[i] = m*(2*p*xc[i] -xc[i]*xc[i])/(p*p)
            dycdxc[i] = 2*m*(p-xc[i])/(p*p)
        else:
            yc[i] = m*(1-2*p+ 2*p*xc[i]-xc[i]*xc[i] )/(1-p)/(1-p)
            dycdxc[i] =  2*m*(p-xc[i])/(1-p)/(1-p)
        
        a0 = 0.2969
        a1 = -0.126
        a2 = -0.3516
        a3 = 0.2843
        a4 = -0.1036
        yt[i] = ((t/0.2)*(a0*xc[i]**0.5 + a1*xc[i] + a2*xc[i]**2 + a3*xc[i]**3+a4*xc[i]**4))
    ##ok 
    theta = np.arctan(dycdxc)
    xu = xc - yt*np.sin(theta)
    yu = yc + yt*np.cos(theta)
    xl = xc+yt*np.sin(theta)
    yl = yc-yt*np.cos(theta)
    
    if normalised:
        # rotation first
        locmin = np.argmin(xu)
        dy = yu[locmin]
        dx = np.min(xu)
        alph = np.arctan(dy/(1-dx))
        if np.min(xl) < dx:
            locmin = np.argmin(xl)
            dy = yu[locmin_x]
            dx = np.min(xl)
            alph = np.arctan(dy/(1-dx))
        ## rotation
        def rotate(cx, cy, x, y, alpha):
            beta =  np.arctan((y-cy)/(cx-x)) - alpha
            d = ((x-cx)**2+(y-cy)**2)**0.5
            xnew = cx - d*np.cos(beta)
            ynew = cy + d*np.sin(beta)
            return xnew, ynew
        
        for i in range(len(xu)):
            x ,y = xu[i], yu[i]
            xu[i], yu[i] = rotate(1.,0,x,y,alph)
        for i in range(len(xl)):
            x ,y = xl[i], yl[i]
            xl[i], yl[i] = rotate(1.,0,x,y,alph)   
            
        minVal = min(np.min(xu),np.min(xl))
        xu = xu - minVal
        xl = xl - minVal
        maxVal = max(np.max(xu), np.max(xl))
        scale = 1./maxVal
        xu *= scale
        xl *= scale
        yu *= scale
        yl *= scale
    return np.stack((xu, yu),axis=1), np.stack((xl, yl),axis=1)
    
        
if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import copy
    myProfile = Aerofoil("myCase",*MPXX_family("9912"))
    oldProfile = copy.deepcopy(myProfile)

    profile= myProfile.getProfile()
    xu,yu, xl, yl = myProfile.getProfile_part(x_min=0.5, x_max=1, nPoints=30)    
    plt.figure()
    plt.plot(profile[:,0], profile[:,1],"r-", label=f"ref {myProfile.Type}")    
    plt.plot(xu,yu, "gv", label="sampling naca")
    plt.plot(xl,yl, "k+", label="sampling naca")
    plt.grid()
        
    xu, yu, xl, yl = myProfile.getProfile_part(x_min=0.1, x_max=1, nPoints=30)
    upperPart = np.stack((xu,yu),axis=1)
    lowerPart = np.stack((xl,yl),axis=1)
    myProfile.update(upperPart,lowerPart)
    profile = myProfile.getProfile()
    xu, yu, xl, yl = myProfile.getProfile_part(x_min=0.7, x_max=1, nPoints=50)
    plt.plot(profile[:,0], profile[:,1],"r.", label=f"profile {myProfile.Type}")
    plt.plot(xu,yu, "mx", label="sampling new")
    plt.plot(xl,yl, "k+", label="sampling new")
    plt.xlabel("x/c")
    plt.ylabel("y/c")
    plt.axis("equal")
    plt.legend()
    
    print(oldProfile.diff(myProfile))
    print(oldProfile.diff_list([(-10,0)]))
    
    plt.figure()
    myProfile = Aerofoil("myCase",*MPXX_family("4412",  normalised=True))
    profile = myProfile.getProfile()
    plt.plot(profile[:,0], profile[:,1],"r-", label=f"ref {myProfile.Type} normalised") 
    myProfile = Aerofoil("myCase",*MPXX_family("4412",  normalised=False))
    profile = myProfile.getProfile()
    plt.plot(profile[:,0], profile[:,1],"g--", label=f"ref {myProfile.Type} orig") 
    plt.grid()
    plt.legend()
    
    plt.show()
    
