# array of particles:
#[x, y, z, radius, mass, pressure, velocity, type]

from numpy import pi, power, array, linalg, sqrt

class Cubic_Spline:
    def __init__(self,r=1.0,h=2.0,dim=3.0,step=1/1e6):
        self.r = r
        self.h = h
        self.step = step
        self.q = abs(r/h)
        self.fac = float(1./power(h,dim))

    def Kernel(self):  
        if self.q >= 0 and self.q < 1:
            return self.fac*3/(2*pi)*(2/3-power(self.q,2)+1/2*power(self.q,3))
        elif self.q >= 1 and self.q <= 2:
            return self.fac/(4*pi)*power(2-self.q,3)
        else:
            return 0
    
    def Gradient(self):
               
        pts = [abs((self.r - self.step)/self.h),abs((self.r + self.step)/self.h)]
        kernel_func = []
        for self.q in pts:
            kernel_func.append(self.Kernel())
        
        return (kernel_func[1]-kernel_func[0])/(2*self.step)
    
    def Laplacian(self):
        
        pts = [abs((self.r - self.step)/(self.h)),abs((self.r)/(self.h)),abs((self.r + self.step)/(self.h))]
        kernel_func = []
        for self.q in pts:
            kernel_func.append(self.Kernel())
            
        return (kernel_func[0] - 2*kernel_func[1] + kernel_func[2])/(self.step**2)
    

class B_Spline:
    def __init__(self,r=1,h=2,dim=3,step=1/1e6):
        self.r = abs(r)
        self.h = h
        self.step = step
        self.q = r/h
        self.fac = 1/power(h,dim)
        
    def Kernel(self):
        self.q = abs(self.q)
        if self.q >= 0 and self.q <= 1:
            return self.fac/pi*(1+(3*power(self.q,3)/4)+(3*power(self.q,2)/2))
        elif self.q > 1 and self.q <= 2:
            return self.fac/pi*(0.25*power(2-self.q,3))
        else:
            return 0
    
    def Gradient(self):
                
        pts = [(self.r - self.step)/self.h, (self.r + self.step)/self.h]
        kernel_func = []
        for self.q in pts:
            kernel_func.append(self.Kernel())
        
        return (kernel_func[1]-kernel_func[0])/(2*self.step)
    
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
 
        if abs(self.r) >= 0 and abs(self.r) <= self.h:
            return 315/(64*pi*self.h**9)*power(self.h**2-abs(self.r)**2,3)
        else:
            return 0

    def Gradient(self):
        
        return -945/(32*pi*self.h**9)*self.r*power(self.h**2-abs(self.r)**2,2)
    
    def Laplacian(self):
            
        return -945/(32*pi*self.h**9)*(self.h**2-abs(self.r)**2)*(3*self.h**2-7*abs(self.r)**2)

class Spiky: #Recommended for pressure
    def __init__(self,r=1,h=2,step=1/1e6):
        self.r = r
        self.h = h
        self.step = step

    def Kernel(self): 
        
        if abs(self.r) >= 0 and abs(self.r) <= self.h:
            return 15/(pi*self.h**6)*power(self.h-(abs(self.r)),3)
        else:
            return 0
        
    def Gradient(self):
               
        return -45/(pi*self.h**6)*self.r/abs(self.r)*power(self.h-abs(self.r),2)
    
    def Laplacian(self):
        
        return -90/(pi*self.h**6*abs(self.r))*(self.h-abs(self.r))*(self.h-2*abs(self.r))

class Viscosity:
    def __init__(self,r=1,h=2,step=1/1e6):
        self.r = r
        self.h = h
        self.step = step
    
    def Kernel(self):
       
        if abs(self.r) >= 0 and abs(self.r) <= self.h:
            return 15/(2*pi*self.h**3)*(-(abs(self.r)**3/(2*self.h**3))+(abs(self.r)**2/self.h**2)+self.h/(2*abs(self.r))-1)
        else:
            return 0
    
    def Gradient(self):
               
        return 15/(2*pi*self.h**3)*self.r*(-3*abs(self.r)/(2*self.h**3)+2/self.h**2-self.h/(2*abs(self.r)**3))
    
    def Laplacian(self):
        
        return 45/(pi*self.h**6)*(self.h-abs(self.r))
    
def Kernel_Correction(neighbors,kernel_name,h):
    A = array([[0,0,0,0],
            [0,0,0,0],
            [0,0,0,0],
            [0,0,0,0]])

    for i in neighbors:
        for r_vector in neighbors[i]['moving_r']:
            r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)
            W = kernel_name(r,h).Kernel()

            if len(r_vector) == 3:

                rx = r_vector[0]
                ry = r_vector[1]
                rz = r_vector[2]
                Aij = array([[1,rx,ry,rz],
                            [rx,rx*rx,rx*ry,rx*rz],
                            [ry,ry*rx,ry*ry,ry*rz],
                            [rz,rz*rx,rz*ry,rz*rz]])
                Aij = Aij * W

                A = [[A[j][k] + Aij[j][k]  for k in range(len(A[0]))] for j in range(len(A))]

    beta = linalg.pinv(A)*[[1],[0],[0],[0]]
    return beta[0]