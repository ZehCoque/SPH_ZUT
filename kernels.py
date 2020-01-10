from numpy import pi, power, sqrt

class Kernel:
    def __init__(self,r,h):
        self.r = r
        self.h = h
    
    def Cubic_Spline(self):
        q = self.r/(2*self.h)
        if q >= 0 and q <= 1/2:
            return 8/pi*(1-6*power(q,2)+6*power(q,3))
        elif q > 1/2 and q <= 1:
            return 16/pi*power(1-q,3)
        else:
            return 0
    
    def B_Spline(self):
        q = self.r/self.h
        if q >= 0 and q <= 1:
            return 1/(pi*power(self.h,3))*(1+(3*power(q,3)/4)+(3*power(q,2)/2))
        elif q > 1 and q <= 2:
            return 1/(pi*power(h,3))*(0.25*power(2-q,3))
        else:
            return 0
    
    def Poly_6(self):
        if r >= 0 and r <= h:
            return 315/(64*pi*h**9)*power(h**2-r**2,3)
        else:
            return 0
    
    def Spiky(self):
        if r >= 0 and r <= h:
            return 15/(pi*h**6)*power(h-r,3)
        else:
            return 0
            
    def Viscosity(self):
        if r >= 0 and r <= h:
            return 15/(2*pi*h**3)*(-r**3/(2*h**3)+r**2/h**2+h/(2*r)-1)
        else:
            return 0
    
coord_1 = [4, 3.5, 0]
coord_2 = [6, 3, 0]

r = sqrt(power(coord_1[0]-coord_2[0],2)+power(coord_1[1]-coord_2[1],2)+power(coord_1[2]-coord_2[2],2))
h = 2.5

W = Kernel(r,h).Viscosity()

print(W)
