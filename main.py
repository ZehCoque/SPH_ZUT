
# count_row = df.shape[0]  # gives number of row count
# count_col = df.shape[1]  # gives number of col count

import pandas as pd
from numpy import sqrt, power

import force_fields
import kernels
import utilities

h = 2.
k = 1.
rho_0 = 0.95
mu = 1
axis_array = ['X','Y','Z'] # auxiliary axis
g = 9.81
radius = 1
lam = 0.4
time = 0
final_time = 10

# Get path to save simulation results
paths = utilities.get_paths("./results/")

#Initialization (arbitrary set of properties)
particles = {'X':[0,0.8,1,0],
             'Y':[0,0,-0.4,0],
             'Z':[0,0,0.25,0.25],
             'X Velocity':[-1.,0.8,1.,3.],
             'Y Velocity':[-0.4,0.5,2.,1.],
             'Z Velocity':[1.,2.,1.5,-2.],
             'Pressure': [0.,0.,0.,0.],
             'Density': [0.,0.,0.,0.],
             'Mass':[1,2,1,1],
             'Type':['Water','Water','Oil','Water']}

particles = pd.DataFrame(particles)

#Saving initial conditions
utilities.save_csv(paths[2],time,particles)
utilities.save_vtk(paths[0],time,particles)

# VTK Group

# Stop when simulation time reaches final time
while time < final_time:

    # Get info from last iteration
    filename = paths[2] + '/sph_' + str(round(time,3)) + '.csv'
    particles = pd.read_csv(filename)
    particles['Neighbors'] = particles.apply(lambda r: [],axis=1) # auxiliary addition
    
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
            r = pd.Series(ri-neighbor_particles[axis])
            neighbor_particles['r' + axis.lower()] = r
        particles.at[i, 'Density'] = force_fields.Density(particles.iloc[i],neighbor_particles,h)
        particles.at[i, 'Pressure'] = k*(power(particles.iloc[i]['Density']/rho_0,7)-1)

    # Third iteration through all particles
    #  -Calculates all force fields (pressure,viscosity and others) for each particle i
    for i in range(0,particles.shape[0]):
        for axis in axis_array:
            particles.at[i, axis + ' Pressure Force'] = force_fields.Pressure(particles.iloc[i],neighbor_particles,h,axis) # A force for each axis
            particles.at[i, axis + ' Viscosity Force'] = mu*force_fields.Viscosity(particles.iloc[i],neighbor_particles,h,axis) # A force for each axis
            
    #print(particles)
    # Calculating total force for each axis
    for axis in axis_array:
        particles[axis + ' Total Force'] = particles[axis + ' Pressure Force'] + particles[axis + ' Viscosity Force']

    # particles['Y Total Force'] = particles['Y Total Force'] - particles['Mass'] * g -< Applying gravity

    # Calculating delta_t
    if particles[['X Velocity','Y Velocity','Z Velocity']].max().max() == 0:
        delta_t = 1e-2
    else:
        delta_t = lam * radius*2 / particles[['X Velocity','Y Velocity','Z Velocity']].max().max()
    
    if delta_t + time > final_time:
        delta_t = final_time - time
    
    # Calculating all new positions for each particle i
    particles_dp1 = pd.DataFrame()
    for axis in axis_array:
        particles_dp1[axis + ' Velocity'] = particles[axis + ' Velocity'] + delta_t * particles[axis + ' Total Force']/particles['Mass']
        particles_dp1[axis] = particles[axis] + delta_t * particles[axis + ' Total Force']/particles['Mass']
        
    particles_dp1 = particles_dp1[['X','Y','Z','X Velocity','Y Velocity','Z Velocity']]
    particles_dp1['Mass'] = particles['Mass']
    particles_dp1['Type'] = particles['Type']
    particles_dp1['Density'] = particles['Density']
    particles_dp1['Pressure'] = particles['Pressure']
    time = time + delta_t    

    #Saving each iteration
    utilities.save_csv(paths[2],time,particles_dp1)
    utilities.save_vtk(paths[0],time,particles_dp1)
    
    #Making vtk group
    utilities.add_to_group(paths[0],time,paths[1])

    #(str(time/final_time * 100) + '%')

utilities.save_group(paths[1])