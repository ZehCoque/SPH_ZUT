
# count_row = df.shape[0]  # gives number of row count
# count_col = df.shape[1]  # gives number of col count

import pandas as pd
from numpy import sqrt, power
import force_fields
import kernels

h = 1.
k = 1.
rho_0 = 0.95
mu = 1
axis_array = ['X','Y','Z'] # auxiliary axis
g = 9.81

#Initialization
particles = {'X':[0,0.8,1,0],
             'Y':[0,0,-0.4,0],
             'Z':[0,0,0.25,0.25],
             'Radius':[1,1,1,1],
             'Mass':[1,1,1,1],
             'Pressure':[.0,.0,.0,.0],
             'Density':[.0,.0,.0,.0],
             'X Velocity':[-1.,0.8,1,2.],
             'Y Velocity':[-0.4,0.5,2,1.],
             'Z Velocity':[1.,2.,1.5,-2.],
             'Type':['Water','Water','Water','Water']}

particles = pd.DataFrame(particles)
particles['Neighbors'] = particles.apply(lambda r: [],axis=1)

# First iteration through all particles
#  -Gets all neighbors for each particle i
for i in range(0,particles.shape[0]):

    # Get the current particle's variables temporarily
    xi = particles.loc[[i],['X']].values[0][0]
    yi = particles.loc[[i],['Y']].values[0][0]
    zi = particles.loc[[i],['Z']].values[0][0]

    # Calculate the distance between points
    r = pd.Series(xi-particles['X'])
    particles['rx'] = r
    r = pd.Series(yi-particles['Y'])
    particles['ry'] = r
    r = pd.Series(zi-particles['Z'])
    particles['rz'] = r
    
    Neighbors = pd.Series((abs(particles['rx']) <= h) & (abs(particles['ry']) <= h) & (abs(particles['rz']) <= h) & (sqrt(particles['rx']**2+particles['ry']**2+particles['rz']**2)!=0))
    
    particles.at[i, 'Neighbors'] = Neighbors[Neighbors].index.values

# Second iteration through all particles
#  -Calculates density and pressure for each particle i
for i in range(0,particles.shape[0]):
    
    indexes = particles.loc[[i],['Neighbors']].values[0][0]
    neighbor_particles = particles.loc[indexes]

    for axis in axis_array:
        # Get the current particle's variables temporarily
        ri = particles.loc[[i],[axis]].values[0][0]

        # Calculate the distance between points (only inside the neighborhood of i)
        r = pd.Series(xi-neighbor_particles[axis])
        neighbor_particles['r' + axis.lower()] = r
    
    particles.at[i, 'Density'] = force_fields.Density(particles.iloc[i],neighbor_particles,h)

    particles.at[i, 'Pressure'] = k*(power(particles.iloc[i]['Density']/rho_0,7)-1)

# Third iteration through all particles
#  -Calculates all force fields (pressure,viscosity and others) for each particle i
for i in range(0,particles.shape[0]):
    for axis in axis_array:
        particles.at[i, axis + ' Pressure Force'] = force_fields.Pressure(particles.iloc[i],neighbor_particles,h,axis) # A force for each axis
        particles.at[i, axis + ' Viscosity Force'] = mu*force_fields.Viscosity(particles.iloc[i],neighbor_particles,h,axis) # A force for each axis

# Calculating total force for each axis
for axis in axis_array:
        particles[axis + ' Total Force'] = particles[axis + ' Pressure Force'] + particles[axis + ' Viscosity Force']

# particles['Y Total Force'] = particles['Y Total Force'] - particles['Mass'] * g -< Applying gravity

# Calculating all new positions for each particle i


    