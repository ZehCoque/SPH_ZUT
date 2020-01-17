
# count_row = df.shape[0]  # gives number of row count
# count_col = df.shape[1]  # gives number of col count

import pandas as pd
from numpy import sqrt, power

import force_fields
import kernels
import utilities
from particle_maker import make_prism as prism

h = 2.
k = 0.5
rho_0 = 1
mu = 1
axis_array = ['X','Y','Z'] # auxiliary axis
g = 9.81
radius = 1
lam = 0.4
time = 0
final_time = 10
iteration = 1

# Get path to save simulation results
paths = utilities.get_paths("./results/")

#Initialization (arbitrary set of properties)
particles = prism(0,0,0,4,4,4,1,10,'Water')

particles = pd.DataFrame(particles)

#Saving initial conditions
utilities.save_csv(paths[2],iteration,particles)
utilities.save_vtk(paths[0],iteration,particles)
utilities.add_to_group(paths[0],iteration,time,paths[1])

# VTK Group

# Stop when simulation time reaches final time
while time < final_time:

    # Get info from last iteration
    filename = paths[2] + '/iter_' + str(iteration)  + '.csv'
    particles = pd.read_csv(filename)
    print(particles)
    particles['Neighbors'] = particles.apply(lambda r: [],axis=1) # auxiliary addition
    
    # First iteration through all particles
    #  -Gets all neighbors for each particle i
    for i in range(0,particles.shape[0]):

        # Get the current particle's variables temporarily
        for axis in axis_array:
            ri = particles.loc[[i],[axis]].values[0][0]

            # Calculate the distance between points (only inside the neighborhood of i)
            r = pd.Series(ri-particles[axis])
            particles['r' + axis.lower()] = r
        
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
        delta_t = 1e-4
    else:
        delta_t = lam * radius*2 / abs(particles[['X Velocity','Y Velocity','Z Velocity']].max().max())
    
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

    iteration = iteration + 1
    #Saving each iteration
    utilities.save_csv(paths[2] + '2',iteration,particles)
    utilities.save_csv(paths[2],iteration,particles_dp1)
    utilities.save_vtk(paths[0],iteration,particles_dp1)
    
    #Making vtk group
    utilities.add_to_group(paths[0],iteration,time,paths[1])

    
    #(str(time/final_time * 100) + '%')

utilities.save_group(paths[1])