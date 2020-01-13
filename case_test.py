
# count_row = df.shape[0]  # gives number of row count
# count_col = df.shape[1]  # gives number of col count

import pandas as pd
from numpy import sqrt
#import force_fields
import kernels

#Initialization
particles = {'X':[0,0.8,-1,2],
             'Y':[0,0,0.4,1],
             'Z':[0,0,0.25,-1],
             'Radius':[1,1,1,1],
             'Mass':[1,1,1,1],
             'Pressure':[1,1,1,1],
             'Density':[1,1,1,1],
             'X Velocity':[-1,0.8,1,2],
             'Y Velocity':[-0.4,0.5,2,1],
             'Z Velocity':[0,0,0,0],
             'Type':['Water','Water','Water','Water']}

h = 1

particles = pd.DataFrame(particles)

for i in range(0,particles.shape[0]):

    # Get the current particle's variables temporarily
    X_tmp = particles.loc[[i],['X']].values[0][0]
    Y_tmp = particles.loc[[i],['Y']].values[0][0]
    Z_tmp = particles.loc[[i],['Z']].values[0][0]
    Type_tmp = particles.loc[[i],['Type']].values[0][0]

    # Calculate the distance between points
    r = pd.Series(sqrt((X_tmp-particles['X'])**2+(Y_tmp-particles['Y'])**2+(Z_tmp-particles['Z'])**2))
    particles['ri-rj'] = r
    
    particles['Neighbor'] = (particles['ri-rj'] > 0) & (particles['ri-rj'] <= h)
    particles['Same Material'] = (particles['Type'] == Type_tmp)
    
    print(particles.loc[[1],['ri-rj']].values[0][0]/h)
    W = kernels.Spiky(particles.loc[[1],['ri-rj']].values[0][0],h).Laplacian()
    
    print(W)

