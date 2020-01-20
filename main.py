
# count_row = df.shape[0]  # gives number of row count
# count_col = df.shape[1]  # gives number of col count

import pandas as pd
from numpy import sqrt, power

import force_fields
import kernels
import utilities
from particle_maker import make_prism as prism

h = 2.1
k = 10
rho_0 = 1000
mu = 1
axis_array = ['X','Y','Z'] # auxiliary axis
g = 9.81
radius = 1
mass = (4/3*radius)**3*rho_0
lam_v = 0.4 # delta_t for velocity parameter
lam_f = 0.2 # delta_t for force parameter
time = 0
final_time = 1
iteration = 1
cube_size = 5

# Get path to save simulation results
paths = utilities.get_paths("./results/")

# Initialization
# Moving Particles (Fluid particles)
# moving = prism(0,2,0,cube_size,cube_size,cube_size,radius,mass,'Water')

# # Boundary particles
# boundary = prism(-4,0,0,-2,10,8,radius,mass,'Boundary')
# boundary = utilities.mergeDict(boundary,prism(32,0,0,2,10,8,radius,mass,'Boundary'))
# boundary = utilities.mergeDict(boundary,prism(-4,-2,0,20,-2,8,radius,mass,'Boundary'))
# boundary = utilities.mergeDict(boundary,prism(-4,-2,-4,20,11,-2,radius,mass,'Boundary'))
# boundary = utilities.mergeDict(boundary,prism(-4,-2,16,20,11,2,radius,mass,'Boundary'))

moving = {'X': [0,0],
        'Y': [4,0],
        'Z': [0,0],
            'X Velocity':[0,0],
            'Y Velocity':[-15,0],
            'Z Velocity':[0,0],
            'Pressure': [0,0],
            'Density': [0,0],
            'Mass':[mass,mass],
            'Type':['Water','Water']}

boundary = {'X': [10],
            'Y': [10],
            'Z': [10],
            'X Velocity':[0],
            'Y Velocity':[0],
            'Z Velocity':[0],
            'Pressure': [0],
            'Density': [0],
            'Mass':[mass],
            'Type':['Boundary']}

moving = pd.DataFrame(moving)
boundary = pd.DataFrame(boundary)

# Saving initial conditions
utilities.save_csv(paths[2],iteration,moving)
utilities.save_moving_vtk(paths[0],iteration,moving)
utilities.save_boundary_vtk(paths[0],boundary)
utilities.add_to_group(paths[0],iteration,time,paths[1])

# Stop when simulation time reaches final time
while time < final_time:

    # Get info from last iteration
    filename = paths[2] + '/iter_' + str(iteration)  + '.csv'
    moving = pd.read_csv(filename)
    moving['Fluid Neighbors'] = moving.apply(lambda r: [],axis=1) # auxiliary addition
    moving['Boundary Neighbors'] = moving.apply(lambda r: [],axis=1) # auxiliary addition
    
    # First iteration through all particles
    #  -Gets all neighbors for each particle i
    for i in range(0,moving.shape[0]):

        for axis in axis_array:
            ri = moving.loc[[i],[axis]].values[0][0]
            r = pd.Series(ri-moving[axis])
            moving['r' + axis.lower()] = r
        
        for axis in axis_array:
            ri = moving.loc[[i],[axis]].values[0][0]
            r = pd.Series(ri-boundary[axis])
            boundary['r' + axis.lower()] = r

        Fluid_Neighbors = pd.Series((abs(moving['rx']) <= h) & (abs(moving['ry']) <= h) & (abs(moving['rz']) <= h))
        Boundary_Neighbors = pd.Series((abs(boundary['rx']) <= h) & (abs(boundary['ry']) <= h) & (abs(boundary['rz']) <= h))
        
        moving.at[i, 'Fluid Neighbors'] = Fluid_Neighbors[Fluid_Neighbors].index.values
        moving.at[i, 'Boundary Neighbors'] = Boundary_Neighbors[Boundary_Neighbors].index.values

    # Second iteration through all particles
    #  -Calculates density and pressure for each particle i
    for i in range(0,moving.shape[0]):
        
        # Moving particles
        indexes = moving.loc[[i],['Fluid Neighbors']].values[0][0]
        neighbor_moving = moving.loc[indexes]
        ri = moving.loc[[i],[axis]].values[0][0]
        for axis in axis_array:
            # Calculate the distance between points (only inside the neighborhood of i)
            r = pd.Series(ri-neighbor_moving[axis])
            neighbor_moving['r' + axis.lower()] = r
        
        # Boundary particles
        indexes = moving.loc[[i],['Boundary Neighbors']].values[0][0]
        neighbor_boundary = boundary.loc[indexes]
        
        for axis in axis_array:           

            # Calculate the distance between points (only inside the neighborhood of i)
            r = pd.Series(ri-neighbor_boundary[axis])
            neighbor_boundary['r' + axis.lower()] = r
    
        moving.at[i, 'Density'] = force_fields.Density(moving.iloc[i],neighbor_moving,h) + force_fields.Density(moving.iloc[i],neighbor_boundary,h)
        moving.at[i, 'Pressure'] = max(k*(power(moving.iloc[i]['Density']/rho_0,7)-1),0) # pressure should always be non-negative

    # Third iteration through all particles
    #  -Calculates all force fields (pressure,viscosity and others) for each particle i
    for i in range(0,moving.shape[0]):

        # Moving particles
        indexes = moving.loc[[i],['Fluid Neighbors']].values[0][0]
        neighbor_moving = moving.loc[indexes]
        ri = moving.loc[[i],[axis]].values[0][0]

        for axis in axis_array:
            # Calculate the distance between points (only inside the neighborhood of i)
            r = pd.Series(ri-neighbor_moving[axis])
            neighbor_moving['r' + axis.lower()] = r

        indexes = moving.loc[[i],['Boundary Neighbors']].values[0][0]
        neighbor_boundary = boundary.loc[indexes]
        
        for axis in axis_array:           
            # Calculate the distance between points (only inside the neighborhood of i)
            r = pd.Series(ri-neighbor_boundary[axis])
            neighbor_boundary['r' + axis.lower()] = r

        for axis in axis_array:
            moving.at[i, axis + ' Pressure Force'] = force_fields.Pressure(moving.iloc[i],neighbor_moving,h,axis)
            moving.at[i, axis + ' Viscosity Force'] = mu*force_fields.Viscosity(moving.iloc[i],neighbor_moving,h,axis)
            moving.at[i, axis + ' Repulsion Force'] = force_fields.Repulsion(moving.iloc[i],neighbor_moving,h,axis) + force_fields.Repulsion(moving.iloc[i],neighbor_boundary,h,axis)
            
    #print(particles)
    # Calculating total force for each axis
    for axis in axis_array:
        moving[axis + ' Total Force'] = moving[axis + ' Pressure Force'] + moving[axis + ' Viscosity Force'] + moving[axis + ' Repulsion Force']

    #moving['Y Total Force'] = moving['Y Total Force'] - moving['Mass'] * g # -< Applying gravity

    # Calculating delta_t
    if abs(moving[['X Velocity','Y Velocity','Z Velocity']]).max().max() == 0:
        delta_t_v = float('inf')
    else:
        delta_t_v = lam_v * radius*2 /abs(moving[['X Velocity','Y Velocity','Z Velocity']]).max().max()
    
    delta_t_f = lam_f * sqrt(radius*2/abs(moving[['X Total Force','Y Total Force','Z Total Force']]).max().max())

    delta_t = min(delta_t_f,delta_t_v)

    if delta_t + time > final_time:
        delta_t = final_time - time
    
    # Calculating all new positions for each particle i
    moving_dp1 = pd.DataFrame()
    for axis in axis_array:
        moving_dp1[axis + ' Velocity'] = moving[axis + ' Velocity'] + delta_t * moving[axis + ' Total Force']/moving['Mass']
        moving_dp1[axis] = moving[axis] + delta_t * moving_dp1[axis + ' Velocity']
        
    moving_dp1 = moving_dp1[['X','Y','Z','X Velocity','Y Velocity','Z Velocity']]
    moving_dp1['Mass'] = moving['Mass']
    moving_dp1['Type'] = moving['Type']
    moving_dp1['Density'] = moving['Density']
    moving_dp1['Pressure'] = moving['Pressure']
    time = time + delta_t    

    iteration = iteration + 1
    #Saving each iteration
    utilities.save_csv(paths[2] + '2',iteration,moving)
    utilities.save_csv(paths[2],iteration,moving_dp1)
    utilities.save_moving_vtk(paths[0],iteration,moving_dp1)
    
    #Making vtk group
    utilities.add_to_group(paths[0],iteration,time,paths[1])

    (str(time/final_time * 100) + '%')

utilities.save_group(paths[1])