
from numpy import sqrt, array, pi, dot, sign, around, isnan
import csv
import force_fields
import kernels
import utilities
import particle_maker
import time as tm
import hashing

# # Initialization

# Making particles

fluid_radius = 0.01 # Radius of each particle
boundary_radius = 0.01
rho_0 = 1000 # Density at rest
mass = 1 # temporary mass of each particle

# Fluid Particles
particles = particle_maker.make_prism2([0,0,0],[0.06,0.08,.06],fluid_radius,mass,'Water')

# Boundary Particles
particles = particle_maker.make_box([-0.1,-0.1,-0.1],[.3,.3,.3],boundary_radius,0,'Stainless Steel',dict_index=len(particles),box=particles)

P = len(particles) # Total number of particles

N = 0
fluid_array = []
boundary_array = []
for i in range(0,P): # Counting fluind particles
    if particles[i]['Type'] == 'Water':
        N += 1
        fluid_array.append(i)
    else:
        boundary_array.append(i)
        
B = P - N # Number of boundary particles

print()
print("#"*40)
print('Initializing with:')
print(str(N) + ' fluid particles and ' + str(B) + ' boundary particles')


#Calculating Kernel Support radius
V = N * fluid_radius**3 * 4/3 * pi
h = (3*V*20/(4*pi*N))**(1/3)
mass = rho_0/kernels.Poly_6(0,h).Kernel()
for i in fluid_array:
    particles[i]['Mass'] = mass

#Hashing
print()
print("#"*40)
print('Hashing')
hash_table_size = hashing.nextPrime(2*P)
for i in range(0,P):
    point = [particles[i]['X'],particles[i]['Y'],particles[i]['Z']]
    hashing.Hashing(h,hash_table_size)._add(point,i)

# hash_table = hashing.Hashing(h,hash_table_size).d

# Getting neighbors:
print("#"*40)
print('Getting neighbors for each particle')

neighborhood = []
for i in range(0,P):
    ri = array([particles[i]['X'],particles[i]['Y'],particles[i]['Z']])

    possible_neighbors = hashing.Hashing(h,hash_table_size).possible_neighbors(ri)

    for j in possible_neighbors:
        rj = array([particles[j]['X'],particles[j]['Y'],particles[j]['Z']])
        r_vector = ri - rj
        r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)

        if  r <= h:
            # t = 1 -> Fluid - Fluid interaction
            # t = 0 -> Fluid - Boundary interaction
            # t = -1 -> Boundary - Boundary interaction
            if particles[i]['Type'] == 'Water' and particles[j]['Type'] == 'Water':
                t = 1
            elif particles[i]['Type'] == 'Water' and particles[j]['Type'] != 'Water':
                t = 0
            else:
                t = -1
            neighborhood.append([i,j,r_vector,r,t])

fluid_fluid = []
boundary_boundary = []
fluid_boundary = []
for num, i in enumerate(neighborhood):
    if i[-1] == 1:
        fluid_fluid.append(num)
    if i[-1] == 0:
        fluid_boundary.append(num)
    if i[-1] == -1:
        boundary_boundary.append(num)

for i in boundary_array:
    particles[i]['psi'] = 0

# Calculating the volume of each boundary particle
print("#"*40)
print('Calculating the volume of each boundary particle')
for i in boundary_boundary:
    r = neighborhood[i][3]
    W = kernels.Poly_6(r,h).Kernel()
    particles[neighborhood[i][0]]['psi'] += rho_0/W

# Used constants
k = 3 # Gas stiffness
mu = 3.5*1.003e-3 # Viscosity coefficient
delta = 0.0728 # Surface tension coefficient
g = array([0,-9.81,0])
lam_v = 0.4 # delta_t for velocity parameter
lam_f = 1 # delta_t for force parameter
gamma = 7 # power of pressure calculation
time = 0 # initial time
final_time = 1 #in seconds
iteration = 1 #First iteration number
start_time = tm.time() # actual clock time of initialization
alpha = 0.01 # First artificial viscosity term
beta_vis = 0 # Second artificial viscosity term
w = 50 # steps for density correction
correction = False # True for correction in the first iteration
delta_friction = 1 # Friction coefficient between 2 interaction surfaces
c_R = 0.01 # coefficient of restitution (between 0 and 1)
beta = []

#Sound velocities (m/s)
c= {'Water':1480,
    'Stainless Steel':5800,
    'Air':330,
    'Polystyrene':2400}

# Get path to save simulation results
paths = utilities.get_paths("./results/")

print("#"*40)
print('Saving initial conditions')
# Saving initial conditions
utilities.save_csv(paths[2],iteration,{key: particles[key] for key in fluid_array})
utilities.save_moving_vtk(paths[0],iteration,{key: particles[key] for key in fluid_array})
try:
    utilities.save_boundary_vtk(paths[0],{key: particles[key] for key in boundary_array})
except:
    pass
utilities.add_to_group(paths[0],iteration,time,paths[1])

print("#"*40)

# Stop when simulation time reaches final time
while time < final_time:

    # Cleaning Force fields
    for i in fluid_array:
        particles[i]['Boundary-Fluid Pressure'] = array([0.,0.,0.])
        particles[i]['Boundary-Fluid Friction'] = array([0.,0.,0.])
        particles[i]['Pressure Force'] = array([0.,0.,0.])
        particles[i]['Viscosity Force'] = array([0.,0.,0.])
        particles[i]['Surface Tension Force'] = array([0.,0.,0.])
        particles[i]['Density'] = 0 

    #  Calculating density for each particle i 
    #  For fluid-fluid interactions: 
    for i in fluid_fluid:
        current_particle = particles[neighborhood[i][0]]
        if correction == True:
            check = neighborhood[i][0]
            count = 0
            neighbors = {}
            j = i
            while neighborhood[i][0] == check:
                neighbors[neighborhood[j][1]] = particles[neighborhood[j][1]]
                count += 1
                j += 1
                try:
                    check = neighborhood[j][0]
                except:
                    break

            ri = array([current_particle['X'],current_particle['Y'],current_particle['Z']])
            beta = kernels.Kernel_Correction(neighbors,ri,h,'Poly_6')
        
        neighbor = particles[neighborhood[i][1]]
        r_vector = neighborhood[i][2]
        r = neighborhood[i][3]
        t = neighborhood[i][4]
        particles[neighborhood[i][0]]['Density'] += force_fields.Density(current_particle,neighbor,r,h,t,r_vector,beta,correction,"Poly_6")
    
    #  For fluid-boundary interactions: 
    for i in fluid_boundary:
        current_particle = particles[neighborhood[i][0]]
        if correction == True:
            check = neighborhood[i][0]
            count = 0
            neighbors = {}
            j = i
            while neighborhood[i][0] == check:
                neighbors[neighborhood[j][1]] = particles[neighborhood[j][1]]
                count += 1
                j += 1
                try:
                    check = neighborhood[j][0]
                except:
                    break

            ri = array([current_particle['X'],current_particle['Y'],current_particle['Z']])
            beta = kernels.Kernel_Correction(neighbors,ri,h,'Poly_6')
        
        neighbor = particles[neighborhood[i][1]]
        r_vector = neighborhood[i][2]
        r = neighborhood[i][3]
        t = neighborhood[i][4]
        particles[neighborhood[i][0]]['Density'] += force_fields.Density(current_particle,neighbor,r,h,t,r_vector,beta,correction,"Poly_6")

    
    # Calculating density for each particle i 
    for i in fluid_array:
        particles[i]['Pressure'] = k*((particles[i]['Density']/rho_0)**gamma-1)
    
    #  -Calculating all force fields (pressure,viscosity and others) for fluid-fluid interactions
    for i in fluid_fluid:
        if neighborhood[i][0] != neighborhood[i][1]:
            current_particle = particles[neighborhood[i][0]]
            neighbor = particles[neighborhood[i][1]]
            r_vector = neighborhood[i][2]
            r = neighborhood[i][3]
            t = neighborhood[i][4]
            particles[neighborhood[i][0]]['Pressure Force'] += force_fields.Pressure(current_particle,neighbor,r,h,t,r_vector,correction,"Spiky",alpha,beta_vis,c)
            particles[neighborhood[i][0]]['Viscosity Force'] += force_fields.Viscosity_Kernel(current_particle,neighbor,r,h,t,correction,"Viscosity",mu)
            particles[neighborhood[i][0]]['Surface Tension Force'] += force_fields.Surface_Tension(current_particle,neighbor,r,h,t,r_vector,correction,"Poly_6",delta)
    
    #-Calculating all force fields (pressure,viscosity and others) for fluid-boundary interactions
    for i in fluid_boundary:
        current_particle = particles[neighborhood[i][0]]
        neighbor = particles[neighborhood[i][1]]
        r_vector = neighborhood[i][2]
        r = neighborhood[i][3]
        t = neighborhood[i][4]
        particles[neighborhood[i][0]]['Boundary-Fluid Pressure'] += force_fields.Boundary_Fluid_Pressure(current_particle,neighbor,r,h,r_vector,correction,'Poly_6')
        particles[neighborhood[i][0]]['Boundary-Fluid Friction'] += force_fields.Boundary_Fluid_Friction(current_particle,neighbor,r,h,r_vector,correction,'Poly_6',delta_friction,c)

    max_force = 0
    # Calculating total force for each axis
    for i in fluid_array:
        particles[i]['Total Force'] = particles[i]['Pressure Force'] + particles[i]['Viscosity Force']\
             + particles[i]['Surface Tension Force'] + rho_0 * g + particles[i]['Boundary-Fluid Pressure']\
                  + particles[i]['Boundary-Fluid Friction']
        for j in range(0,3):
            if abs(particles[i]['Total Force'][j]) > abs(max_force):
                max_force = abs(particles[i]['Total Force'][j])
                      
    max_vel = 0
    for i in fluid_array:
        for j in ['X Velocity','Y Velocity','Z Velocity']:
            if abs(particles[i][j]) > max_vel:
                max_vel = abs(particles[i][j])

    # Calculating delta_t
    if max_vel == 0 and max_force == 0:
        delta_t = 1e-2
    else:
        try:
            delta_t_v = lam_v * fluid_radius*2 /max_vel
        except:
            delta_t_v = float('inf')
        
        try:
            delta_t_f = lam_f * sqrt(fluid_radius*2/max_force)
        except:
            delta_t_f = 1e-2

        delta_t = min(delta_t_f,delta_t_v)

    if delta_t + time > final_time:
        delta_t = final_time - time
    
    # Calculating all new positions for each particle i
    for i in fluid_array:
        count = 0
        for j in ['X','Y','Z']:
            particles[i][j + " Velocity"] = particles[i][j + " Velocity"] + delta_t * particles[i]['Total Force'][count]/particles[i]['Density']
            particles[i][j] = particles[i][j] + delta_t * particles[i][j + ' Velocity']
            count = count + 1

    # New neighborhood calculation
    # Hashing
    hashing.Hashing(h,hash_table_size,d= {})
    for i in range(0,P):
        point = [particles[i]['X'],particles[i]['Y'],particles[i]['Z']]
        hashing.Hashing(h,hash_table_size)._add(point,i)

    # neighborhood search
    neighborhood = []
    colliders = []
    for i in fluid_array:
        ri = array([particles[i]['X'],particles[i]['Y'],particles[i]['Z']])

        possible_neighbors = hashing.Hashing(h,hash_table_size).possible_neighbors(ri)

        for j in possible_neighbors:
            rj = array([particles[j]['X'],particles[j]['Y'],particles[j]['Z']])
            r_vector = ri - rj
            r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)

            if  r <= h:
                # t = 1 -> Fluid - Fluid interaction
                # t = 0 -> Fluid - Boundary interaction
                # t = -1 -> Boundary - Boundary interaction
                if particles[i]['Type'] == 'Water' and particles[j]['Type'] == 'Water':
                    t = 1
                elif particles[i]['Type'] == 'Water' and particles[j]['Type'] != 'Water':
                    t = 0
                else:
                    t = -1
                neighborhood.append([i,j,r_vector,r,t])
                
                # Collision detection

                if r < (fluid_radius + boundary_radius):
                    colliders.append([i,j])
    
    fluid_fluid = []
    fluid_boundary = []
    for num, i in enumerate(neighborhood):
        if i[-1] == 1:
            fluid_fluid.append(num)
        elif i[-1] == 0:
            fluid_boundary.append(num)

    # # Collision Response
    # try:
    #     current_particle = colliders[0][0]
    # except:
    #     pass
  
    # for i in range(0,len(colliders)):

    #     if current_particle == colliders[i][0]:
    #         fluid_collider = moving[colliders[i][0]]
    #         boundary_collider = boundary[colliders[i][1]]

    #         r_index = fluid_collider['Boundary Neighbors'].index(colliders[i][1])
    #         string = ['X Velocity','Y Velocity','Z Velocity']
    #         d = abs(sqrt(fluid_collider['boundary_r'][r_index][0]**2+fluid_collider['boundary_r'][r_index][1]**2+fluid_collider['boundary_r'][r_index][2]**2) - (fluid_radius + boundary_radius))
    #         current_particle = colliders[i][0]
    #         new_vel = array([0.,0.,0.])
    #         mag_old_vel = sqrt(fluid_collider[string[0]]**2+fluid_collider[string[1]]**2+fluid_collider[string[2]]**2)
    #         direction = sign(around(fluid_collider['boundary_r'][r_index], decimals = 10))
    #         normal = -1*direction
            
    #         old_vel = []
    #         for v_axis in string:
    #             old_vel.append(fluid_collider[v_axis])
            
    #         new_vel = old_vel - (1+c_R*d/(delta_t*mag_old_vel))*dot(old_vel,normal)*normal

    #         count = 0
    #         for v_axis in string:
    #             moving[colliders[i][0]][v_axis] = new_vel[count] * lam_v * fluid_radius*2 /max(new_vel)
    #             moving[colliders[i][0]][v_axis[0]] = moving[colliders[i][0]][v_axis[0]] + lam_v * fluid_radius*2 /max(new_vel) * new_vel[count]
    #             count += 1
             
    time = time + delta_t    

    iteration = iteration + 1
    correction = False
    # if iteration % w == 0:
    #     correction = True
    
    #Saving iterations
    utilities.save_csv(paths[2],iteration,{key: particles[key] for key in fluid_array})
    utilities.save_moving_vtk(paths[0],iteration,{key: particles[key] for key in fluid_array})
    
    #Making vtk group
    utilities.add_to_group(paths[0],iteration,time,paths[1])

    utilities.info(time,final_time,start_time,tm.time(),delta_t,iteration)

utilities.save_group(paths[1])