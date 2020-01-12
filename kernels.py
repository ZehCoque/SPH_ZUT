# array of particles:
#[x, y, z, radius, mass, pressure, velocity, type]

from numpy import pi, power, sqrt, linspace, insert
import pandas as pd

class Cubic_Spline:
    def __init__(self,posi=[0,0,0],posj=[0,0,0],h=2,grad_dim='x',step=1/1e6):
        self.posi = posi
        self.posj = posj
        self.h = h
        self.step = step
        self.grad_dim = grad_dim

    def Gradient(self):
        if self.grad_dim == 'x':
            dim = 0
        elif self.grad_dim == 'y':
            dim = 1
        elif self.grad_dim == 'z':
            dim = 2
        
        r = self.posi[0][dim] - self.posj[0][dim]

        pts = [(r - self.step)/(self.h),(r + self.step)/(self.h)]
        kernel_func = []
        for q in pts:
            if q >= 0 and q <= 1/2:
                kernel_func.append(8/pi*(1-6*power(q,2)+6*power(q,3)))
            elif q > 1/2 and q <= 1:
                kernel_func.append(16/pi*power(1-q,3))
            else:
                kernel_func.append(0)
        
        return (kernel_func[1]-kernel_func[0])/(2*self.step)

    def Kernel(self):
        xi = self.posi[0]
        xj = self.posj[0]
        yi = self.posi[1]
        yj = self.posj[1]
        zi = self.posi[2]
        zj = self.posj[2]
        r = sqrt(power(xi-xj,2)+power(yi-yj,2)+power(zi-zj,2))
        q = r/(self.h*2)
        if q >= 0 and q <= 1/2:
            return 8/pi*(1-6*power(q,2)+6*power(q,3))
        elif q > 1/2 and q <= 1:
            return 16/pi*power(1-q,3)
        else:
            return 0

class B_Spline:
    def __init__(self,r=1.,h=2.):
        self.r = r
        self.h = h
    
    def Kernel(self):
        q = self.r/self.h
        if q >= 0 and q <= 1:
            return 1/(pi*power(self.h,3))*(1+(3*power(q,3)/4)+(3*power(q,2)/2))
        elif q > 1 and q <= 2:
            return 1/(pi*power(self.h,3))*(0.25*power(2-q,3))
        else:
            return 0

class Poly_6:
    def __init__(self,r=1,h=2,step=1/1e6):
        self.r = r
        self.h = h
        self.step = step

    def Gradient(self):
        pts = [self.r - self.step, self.r + self.step]

        if pts[0] < 0:
            pts[0]=0
        if pts[1] > self.h:
            pts[1] = self.h

        kernel_func = []
        for r in pts:
            if r >= 0 and r <= self.h:
                kernel_func.append(315/(64*pi*self.h**9)*power(self.h**2-r**2,3))
        
        return abs((kernel_func[1]-kernel_func[0])/(2*self.step))

    def Kernel(self):
        if self.r >= 0 and self.r <= self.h:
            return 315/(64*pi*self.h**9)*power(self.h**2-self.r**2,3)
        else:
            return 0

class Spiky:
    def __init__(self,r=1,h=2):
        self.r = r
        self.h = h

    def Kernel(self): #Recommended for pressure
        if self.r >= 0 and self.r <= self.h:
            return 15/(pi*self.h**6)*power(self.h-self.r,3)
        else:
            return 0

class Viscosity:
    def __init__(self,r=1,h=2):
        self.r = r
        self.h = h
            
    def Kernel(self):
        if self.r >= 0 and self.r <= self.h:
            return 15/(2*pi*self.h**3)*(-self.r**3/(2*self.h**3)+self.r**2/self.h**2+self.h/(2*self.r)-1)
        else:
            return 0


# coord_1 = [4, 3.5, 0]
# coord_2 = [6, 3, 0]

# r = sqrt(power(coord_1[0]-coord_2[0],2)+power(coord_1[1]-coord_2[1],2)+power(coord_1[2]-coord_2[2],2))
# r = 0.5
# h = 1

# Gradient = Poly_6(r,h).Gradient()

# print(Gradient)