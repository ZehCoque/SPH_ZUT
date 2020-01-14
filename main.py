
# count_row = df.shape[0]  # gives number of row count
# count_col = df.shape[1]  # gives number of col count

import pandas as pd
from numpy import sqrt
#import force_fields
import kernels

#Initialization
particles = {'X':[0,0.8,1,0],
             'Y':[0,0,-0.4,0],
             'Z':[0,0,0.25,0.25],
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
particles['Neighbors'] = particles.apply(lambda r: [],axis=1)

# First iteration through all particles
#  -Gets all neighbors for each particle i
for i in range(0,particles.shape[0]):

    # Get the current particle's variables temporarily
    xi = particles.loc[[i],['X']].values[0][0]
    yi = particles.loc[[i],['Y']].values[0][0]
    zi = particles.loc[[i],['Z']].values[0][0]
    Type_i = particles.loc[[i],['Type']].values[0][0]

    # Calculate the distance between points
    r = pd.Series(xi-particles['X'])
    particles['rx'] = r
    r = pd.Series(yi-particles['Y'])
    particles['ry'] = r
    r = pd.Series(zi-particles['Z'])
    particles['rz'] = r
    
    Neighbors = pd.Series((abs(particles['rx']) <= h) & (abs(particles['ry']) <= h) & (abs(particles['rz']) <= h) & (sqrt(particles['rx']**2+particles['ry']**2+particles['rz']**2)!=0))
    
    particles.at[i, 'Neighbors'] = Neighbors[Neighbors].index.values

    
    print(particles)
    
    #particles['Same Material'] = (particles['Type'] == Type_i)
    # print(particles.loc[[1],['ri-rj']].values[0][0]/h)
    # W = kernels.Spiky(particles.loc[[1],['ri-rj']].values[0][0],h).Laplacian()
    
    #print(W)

