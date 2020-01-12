
# from numpy import mgrid, gradient

# x,y,z = mgrid[-100:101:25., -100:101:25., -100:101:25.]

# V = 2*x**2 + 3*y**2 - 4*z # just a random function for the potential

# Ex,Ey,Ez = gradient(V)

# print(V)

import pandas as pd
import kernels

particles = {'X':[4,6],
             'Y':[6,3],
             'Z':[0,0],
             'Radius':[1,1],
             'Mass':[1,1],
             'Pressure':[1,1],
             'X Velocity':[-1,0.8],
             'Y Velocity':[-0.4,0.5],
             'Z Velocity':[0,0],
             'Type':['Water','Water']}

particles = pd.DataFrame(particles)

posi = particles.loc[[0],['X','Y','Z']].to_numpy()
posj = particles.loc[[1],['X','Y','Z']].to_numpy()
h = 2

Gradient = kernels.Cubic_Spline(posi,posj,h,'x').Gradient()

print(Gradient)