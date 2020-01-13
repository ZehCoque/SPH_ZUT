# array of particles:
#[x, y, z, radius, mass, pressure, velocity, type]

from numpy import pi, power, sqrt, linspace, insert
import pandas as pd

class Cubic_Spline:
    def __init__(self,r=1,h=2,step=1/1e6):
        self.r = r
        self.h = h
        self.step = step

    def Kernel(self):
                
        if self.q >= 0 and self.q <= 1/2:
            return 8/pi*(1-6*power(self.q,2)+6*power(self.q,3))
        elif self.q > 1/2 and self.q <= 1:
            return 16/pi*power(1-self.q,3)
        else:
            return 0
    
    def Gradient(self):
               
        pts = [(self.r - self.step)/(2*self.h),(self.r + self.step)/(2*self.h)]
        kernel_func = []
        for self.q in pts:
            kernel_func.append(self.Kernel())
        
        return -(kernel_func[1]-kernel_func[0])/(2*self.step)
    
    def Laplacian(self):
        
        pts = [(self.r - self.step)/(2*self.h),(self.r)/(2*self.h),(self.r + self.step)/(2*self.h)]
        kernel_func = []
        for self.q in pts:
            kernel_func.append(self.Kernel())
            
        return 2*kernel_func[0]/(2*self.step**2) - 2*kernel_func[1]/(self.step**2) + 2*kernel_func[2]/(2*self.step**2)
    

class B_Spline:
    def __init__(self,r=1,h=2,step=1/1e6):
        self.r = r
        self.h = h
        self.step = step
        
    def Kernel(self):
        
        if q >= 0 and q <= 1:
            return 1/(pi*power(self.h,3))*(1+(3*power(q,3)/4)+(3*power(q,2)/2))
        elif q > 1 and q <= 2:
            return 1/(pi*power(self.h,3))*(0.25*power(2-q,3))
        else:
            return 0
    
    def Gradient(self):
               
        pts = [(self.r - self.step)/self.h, (self.r + self.step)/self.h]
        kernel_func = []
        for self.q in pts:
            kernel_func.append(self.Kernel())
        
        return -(kernel_func[1]-kernel_func[0])/(2*self.step)
    
    def Laplacian(self):
        
        pts = [(self.r - self.step)/self.h, (self.r)/self.h, (self.r + self.step)/self.h]
        kernel_func = []
        for self.q in pts:
            kernel_func.append(self.Kernel())
            
        return (kernel_func[0] - 2*kernel_func[1] + kernel_func[2])/(self.step**2)  

class Poly_6:
    def __init__(self,r=1,h=2,step=1/1e6):
        self.r = r
        self.h = h
        self.step = step
        
    def Kernel(self):
 
        if self.r >= 0 and self.r <= self.h:
            return 315/(64*pi*self.h**9)*power(self.h**2-self.r**2,3)
        else:
            return 0

    def Gradient(self):
               
        pts = [(self.r - self.step), (self.r + self.step)]
        kernel_func = []
        for self.r in pts:
            kernel_func.append(self.Kernel())
        
        return -(kernel_func[1]-kernel_func[0])/(2*self.step)
    
    def Laplacian(self):
        
        pts = [(self.r - self.step),(self.r),(self.r + self.step)]
        kernel_func = []
        for self.r in pts:
            kernel_func.append(self.Kernel())
            
        return (kernel_func[0] - 2*kernel_func[1] + kernel_func[2])/(self.step**2)

class Spiky: #Recommended for pressure
    def __init__(self,r=1,h=2,step=1/1e6):
        self.r = r
        self.h = h
        self.step = step

    def Kernel(self): 
        
        if self.r >= 0 and self.r <= self.h:
            return 15/(pi*self.h**6)*power(self.h-self.r,3)
        else:
            return 0
        
    def Gradient(self):
               
        pts = [(self.r - self.step), (self.r + self.step)]
        kernel_func = []
        for self.r in pts:
            kernel_func.append(self.Kernel())
        
        return -(kernel_func[1]-kernel_func[0])/(2*self.step)
    
    def Laplacian(self):
        
        pts = [(self.r + self.step), (self.r), (self.r - self.step)]
        kernel_func = []
        for self.r in pts:
            kernel_func.append(self.Kernel())
            
        return (kernel_func[0] - 2*kernel_func[1] + kernel_func[2])/(self.step**2)

class Viscosity:
    def __init__(self,r=1,h=2,step=1/1e6):
        self.r = r
        self.h = h
        self.step = step
    
    def Kernel(self):
       
        if self.r >= 0 and self.r <= self.h:
            return 15/(2*pi*self.h**3)*(-(self.r**3/(2*self.h**3))+(self.r**2/self.h**2)+self.h/(2*self.r)-1)
        else:
            return 0
    
    def Gradient(self):
               
        pts = [(self.r - self.step), (self.r + self.step)]
        kernel_func = []
        for self.r in pts:
            kernel_func.append(self.Kernel())
        
        return -(kernel_func[1]-kernel_func[0])/(2*self.step)
    
    def Laplacian(self):
        
        pts = [(self.r - self.step), (self.r), (self.r + self.step)]
        kernel_func = []
        for self.r in pts:
            kernel_func.append(self.Kernel())
            
        return (kernel_func[0] - 2*kernel_func[1] + kernel_func[2])/(self.step**2)  