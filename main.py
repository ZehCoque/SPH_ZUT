
# count_row = df.shape[0]  # gives number of row count
# count_col = df.shape[1]  # gives number of col count
from numpy import sqrt, power, array
import csv
import force_fields
import kernels
import utilities
from particle_maker import make_prism as prism

h = 0.2
k = 2000
rho_0 = 1000
mu = 3.5
g = array([0,-9.81,0])
radius = h/2
mass = (4/3*radius)**3*rho_0
lam_v = 0.4 # delta_t for velocity parameter
lam_f = 0.2 # delta_t for force parameter
time = 0
final_time = 1
iteration = 1
cube_size = 3

# Get path to save simulation results
paths = utilities.get_paths("./results/")

# # Initialization
# # Boundary particles
# boundary = prism(-4,0,0,-2,10,8,radius,mass,'Boundary')
# boundary = prism(32,0,0,2,10,8,radius,mass,'Boundary',count=len(boundary),prism=boundary)
# boundary = prism(-4,-2,0,20,-2,8,radius,mass,'Boundary',count=len(boundary),prism=boundary)
# boundary = prism(-4,-2,-4,20,11,-2,radius,mass,'Boundary',count=len(boundary),prism=boundary)
# boundary = prism(-4,-2,16,20,11,2,radius,mass,'Boundary',count=len(boundary),prism=boundary)

# # Moving Particles (Fluid particles)
# moving = prism(0,8,0,cube_size,cube_size,cube_size,radius,mass,'Water',count=0,prism={})

moving = {}
boundary = {}

moving[0] = {'X': 0.,'Y': -0.05,'Z': 0.,'X Velocity':0.,'Y Velocity':0.,
            'Z Velocity':0.,'Pressure': 0.,'Density': 0.,'Mass':mass,'Type':'Water'}
moving[1] = {'X': 0.,'Y': 0.05,'Z': 0.,'X Velocity':0.,'Y Velocity':0.,
            'Z Velocity':0.,'Pressure': 0.,'Density': 0.,'Mass':mass,'Type':'Water'}
moving[2] = {'X': 0.05,'Y': 0,'Z': 0.,'X Velocity':0.,'Y Velocity':0.,
            'Z Velocity':0.,'Pressure': 0.,'Density': 0.,'Mass':mass,'Type':'Water'}
moving[3] = {'X': -0.05,'Y': 0,'Z': 0.,'X Velocity':0.,'Y Velocity':0.,
            'Z Velocity':0.,'Pressure': 0.,'Density': 0.,'Mass':mass,'Type':'Water'}

boundary[0] = {'X': 0.,'Y': 0.,'Z': 0.,'X Velocity':0.,'Y Velocity':0.,
            'Z Velocity':0.,'Pressure': 0.,'Density': 0.,'Mass':mass,'Type':'Water'}

# Saving initial conditions
utilities.save_csv(paths[2],iteration,moving)
utilities.save_moving_vtk(paths[0],iteration,moving)
utilities.save_boundary_vtk(paths[0],boundary)
utilities.add_to_group(paths[0],iteration,time,paths[1])

# Stop when simulation time reaches final time
while time < final_time:
    
    # First iteration through all particles
    #  -Gets moving neighbors for each moving particle i
    for i in range(0,len(moving)):
        moving[i]['Fluid Neighbors'] = []
        moving[i]['moving_r'] = []
        ri_x = moving[i]['X']
        ri_y = moving[i]['Y']
        ri_z = moving[i]['Z']

        for j in range(0,len(moving)):
                rx = ri_x-moving[j]['X']
                ry = ri_y-moving[j]['Y']
                rz = ri_z-moving[j]['Z']

                if abs(rx) <= h and abs(ry) <= h and abs(rz) <= h:
                    moving[i]['moving_r'].append([rx,ry,rz])
                    moving[i]['Fluid Neighbors'].append(j)

    #  -Gets boundary neighbors for each moving particle i
    for i in range(0,len(moving)):
        moving[i]['Boundary Neighbors'] = []
        moving[i]['boundary_r'] = []
        ri_x = moving[i]['X']
        ri_y = moving[i]['Y']
        ri_z = moving[i]['Z']
        for j in range(0,len(boundary)):
                rx = ri_x-boundary[j]['X']
                ry = ri_y-boundary[j]['Y']
                rz = ri_z-boundary[j]['Z']

                if abs(rx) <= h and abs(ry) <= h and abs(rz) <= h:
                    moving[i]['boundary_r'].append([rx,ry,rz])
                    moving[i]['Boundary Neighbors'].append(j)

    # Second iteration through all particles
    #  -Calculates density and pressure for each particle i
    for i in range(0,len(moving)):
        moving[i]['Density'] = 0.
        moving[i]['Density'] = force_fields.Density(moving,boundary,i,h)
        moving[i]['Pressure'] = 0.
        moving[i]['Pressure'] = k*(power(moving[i]['Density']/rho_0,7)-1) # pressure should always be non-negative

    # Third iteration through all particles
    #  -Calculates all force fields (pressure,viscosity and others) for each particle i
    for i in range(0,len(moving)):
        moving[i]['Pressure Force'] = force_fields.Pressure(moving,boundary,i,h)
        moving[i]['Viscosity Force'] = force_fields.Viscosity(moving,i,h,mu)
        #moving[i]['Repulsion Force'] = force_fields.Repulsion(moving.iloc[i],neighbor_moving,h,axis) + force_fields.Repulsion(moving.iloc[i],neighbor_boundary,h,axis)

    max_force = 0
    # Calculating total force for each axis
    for i in range(0,len(moving)):
            moving[i]['Total Force'] = moving[i]['Pressure Force'] + moving[i]['Viscosity Force'] #+ moving[i]['Mass'] * g/100 #+ moving[axis + ' Repulsion Force']
            for j in range(0,3):
                if moving[i]['Total Force'][j] > abs(max_force):
                    max_force = abs(moving[i]['Total Force'][j])
                    
    max_vel = 0
    for i in range(0,len(moving)):
        for j in ['X Velocity','Y Velocity','Z Velocity']:
            if abs(moving[i][j]) > max_vel:
                max_vel = abs(moving[i][j])

    # Calculating delta_t
    if max_vel == 0 and max_force == 0:
        delta_t = 1e-2
    else:
        try:
            delta_t_v = lam_v * radius*2 /max_vel
        except:
            delta_t_v = float('inf')
        
        try:
            delta_t_f = lam_f * sqrt(radius*2/max_force)
        except:
            delta_t_f = 1e-2
            
        delta_t = min(delta_t_f,delta_t_v)

    if delta_t + time > final_time:
        delta_t = final_time - time
    
    # Calculating all new positions for each particle i
    for i in range(0,len(moving)):
        count = 0
        for j in ['X','Y','Z']:
            moving[i][j + " Velocity"] = moving[i][j + " Velocity"] + delta_t * moving[i]['Total Force'][count]/moving[i]['Mass']
            moving[i][j] = moving[i][j] + delta_t * moving[i][j + ' Velocity']
            count = count + 1

    time = time + delta_t    

    iteration = iteration + 1
    #Saving each iteration
    utilities.save_csv(paths[2],iteration,moving)
    utilities.save_moving_vtk(paths[0],iteration,moving)
    
    #Making vtk group
    utilities.add_to_group(paths[0],iteration,time,paths[1])

    print(str(round(time/final_time * 100,3)) + '%')

utilities.save_group(paths[1])