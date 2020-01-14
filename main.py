
# count_row = df.shape[0]  # gives number of row count
# count_col = df.shape[1]  # gives number of col count

import pandas as pd
from numpy import sqrt
import force_fields
import kernels

#Initialization
particles = {'X':[0,0.8,1,0],
             'Y':[0,0,-0.4,0],
             'Z':[0,0,0.25,0.25],
             'Radius':[1,1,1,1],
             'Mass':[1,1,1,1],
             'Pressure':[1000,1000,1000,1000],
             'X Density':[1000,1000,1000,1000],
             'Y Density':[1000,1000,1000,1000],
             'Z Density':[1000,1000,1000,1000],
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

    # print(particles)
    
for i in range(0,particles.shape[0]):
    
    indexes = particles.loc[[i],['Neighbors']].values[0][0]
    neighbor_particles = particles.loc[indexes]

    density_xi = particles.loc[[i],['X Density']].values[0][0]
    density_yi = particles.loc[[i],['Y Density']].values[0][0]
    density_zi = particles.loc[[i],['Z Density']].values[0][0]
    r_xi = particles.loc[[i],['rx']].values[0][0]
    r_yi = particles.loc[[i],['ry']].values[0][0]
    r_zi = particles.loc[[i],['rz']].values[0][0]
    
    Grad_W_x = kernels.Cubic_Spline(r=r_xi,h=h)
    Grad_W_y = kernels.Cubic_Spline(r=r_yi,h=h)
    Grad_W_z = kernels.Cubic_Spline(r=r_zi,h=h)

    particles['X Density'] = force_fields.Density(neighbor_particles['Mass'],density_xi,neighbor_particles['X Density'],Grad_W_x)
    particles['Y Density'] = force_fields.Density(neighbor_particles['Mass'],density_yi,neighbor_particles['X Density'],Grad_W_y)
    particles['Z Density'] = force_fields.Density(neighbor_particles['Mass'],density_zi,neighbor_particles['X Density'],Grad_W_z)

    print(particles)

    #particles['Same Material'] = (particles['Type'] == Type_i)
    # print(particles.loc[[1],['ri-rj']].values[0][0]/h)
    # W = kernels.Spiky(particles.loc[[1],['ri-rj']].values[0][0],h).Laplacian()
    
    #print(W)

