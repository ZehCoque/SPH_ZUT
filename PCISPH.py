
from numpy import sqrt,dot, array, pi, dot, sign, around, floor, sign
import force_fields
import kernels
import utilities
import particle_maker
import time as tm
import hashing

# # Initialization

fluid_radius = 0.01 # Radius of each particle
boundary_radius = 0.01
rho_0 = 1000 # Density at rest

# Continuing from last simulation
continue_from_last = False
if continue_from_last == True:
    pull_last = utilities.continue_last_sim()
    particles = pull_last[0]
    particles = particle_maker.make_box([-0.1,-0.1,-0.1],[.3,.3,.3],boundary_radius,0,0,'Stainless Steel',dict_index=len(particles),box=particles)
    iteration = pull_last[1] # Last iteration number
    with open(pull_last[2] + "/time.txt", "r") as f:
        lines = f.read().splitlines()
        time = float(lines[-1])

else:
    # Making particles
    mass = 1 # temporary mass of each particle
    iteration = 1 #First iteration number
    time = 0 # initial time

    # Fluid Particles
    particles = particle_maker.make_prism2([-0.04,-0.04,-0.04],[0.2,0.24,0.2],fluid_radius,mass,rho_0,'Water')

    # Boundary Particles
    particles = particle_maker.make_box([-0.1,-0.1,-0.1],[.3,.3,.3],boundary_radius,0,0,'Stainless Steel',dict_index=len(particles),box=particles)

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
x = 45
h = fluid_radius*4

#Calculating the mass of each particle
if continue_from_last == False:
    mass = 4/3*pi*fluid_radius**3*rho_0
    for i in fluid_array:
        particles[i]['Mass'] = mass

#Hashing
hash_table_size = hashing.nextPrime(P*P)
for i in range(0,P):
    point = [particles[i]['X'],particles[i]['Y'],particles[i]['Z']]
    hashing.Hashing(h,hash_table_size)._add(point,i)

# hash_table = hashing.Hashing(h,hash_table_size).d

# Getting neighbors:
print("#"*40)
print('Getting neighbors for each particle')

boundary_boundary = []
combinations = []
for i in boundary_array:
    ri = array([particles[i]['X'],particles[i]['Y'],particles[i]['Z']])

    possible_neighbors = hashing.Hashing(h,hash_table_size).possible_neighbors(ri)

    for j in possible_neighbors:
        combinations.append([i,j])
        rj = array([particles[j]['X'],particles[j]['Y'],particles[j]['Z']])
        r_vector = ri - rj
        r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)

        if  r < h:
            t = -1
            boundary_boundary.append([i,j,r_vector,r,t])

for i in boundary_array:
    particles[i]['psi'] = 0

# Calculating the volume of each boundary particle
print("#"*40)
print('Calculating the volume of each boundary particle')
for i in boundary_boundary:
    r = i[3]
    W = kernels.Poly_6(r,h).Kernel()
    particles[i[0]]['psi'] += rho_0/W

# Used constants
mu = 3.5*1.003e-3 # Viscosity coefficient
delta = 0.0728 # Surface tension coefficient
g = array([0,-9.81,0])
gamma = 7 # power of pressure calculation
alpha = 0.01 # First artificial viscosity term
beta_vis = 0 # Second artificial viscosity term
delta_friction = 0.01  # Friction coefficient between 2 interaction surfaces
w = 10 # Iterations for density correction
correction = True # True for correction in the first iteration
c_R = 0.01 # coefficient of restitution (between 0 and 1)
final_time = 1 #in seconds

start_time = tm.time() # actual clock time of initialization
beta = [] # cleaning kernel correction beta vector

#PCISPH variables
density_fluctuation_permited = rho_0/100 # 1%
maximum_volume_fluctuation = 10*density_fluctuation_permited
delta_t = 1e-3 # initial delta_t 
# delta_t = 0.25*h/max_vel if particles have initial velocity

# Calculating pressure delta value
ppd = int(floor((h-fluid_radius)/(2*fluid_radius))+1) #particles per direction
coord_init = -2*fluid_radius*ppd
project_particle_array = []
num = ppd*2

particles_array = [float(coord_init)]
for i in range(0,abs(num)):
    particles_array.append(float(particles_array[-1]+2*fluid_radius))

for i in particles_array:
    for j in particles_array:
        for k in particles_array:
            project_particle_array.append([i,j,k])

Grad_W = array([0.,0.,0.])
dot_Grad_W = 0
for i in array(project_particle_array):
    r_vector = -i
    r = sqrt(i[0]**2 + i[1]**2 + i[2]**2 )
    direction = sign(around(r_vector,10))
    Grad_W +=  direction * kernels.Poly_6(r,h).Gradient()
    dot_Grad_W += dot(Grad_W,Grad_W) 

pressure_delta = -(dot(Grad_W,Grad_W)-dot_Grad_W)**-1

#Sound velocities (m/s)
c= {'Water':1480,
    'Stainless Steel':5800,
    'Air':330,
    'Polystyrene':2400}


if continue_from_last==False:
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
else:
    paths = [pull_last[2] + "/vtk", pull_last[3],pull_last[2] + "/csv"]

print("#"*40)
print("Initiating main loop")

# Stop when simulation time reaches final time
while time < final_time:

    # New neighborhood calculation
    # Hashing
    for i in fluid_array:
        point = [particles[i]['X'],particles[i]['Y'],particles[i]['Z']]
        hashing.Hashing(h,hash_table_size)._add(point,i)

    # neighborhood search
    neighborhood = []
    for i in fluid_array:
        ri = array([particles[i]['X'],particles[i]['Y'],particles[i]['Z']])

        possible_neighbors = hashing.Hashing(h,hash_table_size).possible_neighbors(ri)

        for j in possible_neighbors:
            rj = array([particles[j]['X'],particles[j]['Y'],particles[j]['Z']])
            r_vector = ri - rj
            r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)

            if  r < h:
                # t = 1 -> Fluid - Fluid interaction
                # t = 0 -> Fluid - Boundary interaction
                # t = -1 -> Boundary - Boundary interaction
                if particles[i]['Type'] == 'Water' and particles[j]['Type'] == 'Water':
                    t = 1
                elif (particles[i]['Type'] == 'Water' and particles[j]['Type'] != 'Water'):
                    t = 0
                neighborhood.append([i,j,r_vector,r,t])

    # assigning arrays with interaction type
    fluid_fluid = []
    fluid_boundary = []
    for num, i in enumerate(neighborhood):
        if i[-1] == 1:
            fluid_fluid.append(num)
        elif i[-1] == 0 and neighborhood[i[0]]['Type'] == 'Water':
            fluid_boundary.append(num)
    
    # Cleaning Force fields and predicted velocites and positions
    for i in fluid_array:
        particles[i]['Boundary-Fluid Pressure'] = array([0.,0.,0.])
        particles[i]['Boundary-Fluid Friction'] = array([0.,0.,0.])
        particles[i]['Pressure Force'] = array([0.,0.,0.])
        particles[i]['Viscosity Force'] = array([0.,0.,0.])
        particles[i]['Surface Tension Force'] = array([0.,0.,0.])
        particles[i]['Pressure'] = 0

    # Computing viscosity and surface tension forces    
    for i in fluid_fluid:
        if neighborhood[i][0] != neighborhood[i][1]:
            current_particle = particles[neighborhood[i][0]]
            neighbor = particles[neighborhood[i][1]]
            r_vector = neighborhood[i][2]
            r = neighborhood[i][3]
            t = neighborhood[i][4]
            particles[neighborhood[i][0]]['Viscosity Force'] += force_fields.Viscosity_Kernel(current_particle,neighbor,r,h,t,"Viscosity",mu)
            particles[neighborhood[i][0]]['Surface Tension Force'] += force_fields.Surface_Tension(current_particle,neighbor,r,h,t,r_vector,"Poly_6",delta)
    
    # Calculating summation of all non-pressure forces
    for i in fluid_array:
        particles[i]['Total Force'] = particles[i]['Viscosity Force'] + particles[i]['Surface Tension Force'] + rho_0 * g

    # Getting max velocity to calculate delta t
    max_vel = 0
    for i in fluid_array:
        for j in ['X Velocity','Y Velocity','Z Velocity']:
            if abs(particles[i][j]) > max_vel:
                max_vel = abs(particles[i][j])


    max_rho_err = 0
    k=0
    while max_rho_err > density_fluctuation_permited or k < 3:
        
        # Predicting velocity and position (t+1)
        for i in fluid_array:
            count = 0
            for j in ['X','Y','Z']:
                particles[i][j + " Pred Velocity"] = particles[i][j + " Velocity"] + delta_t * particles[i]['Total Force'][count]/particles[i]['Density']
                particles[i][j + " Pred"] = particles[i][j] + delta_t * particles[i][j + ' Velocity']
                count = count + 1
        
    # New neighborhood calculation with predicted positions
    # Hashing
    for i in fluid_array:
        point = [particles[i]['X Pred'],particles[i]['Y Pred'],particles[i]['Z Pred']]
        hashing.Hashing(h,hash_table_size)._add(point,i)

    # neighborhood search
    neighborhood = []
    for i in fluid_array:
        ri = array([particles[i]['X Pred'],particles[i]['Y Pred'],particles[i]['Z Pred']])

        possible_neighbors = hashing.Hashing(h,hash_table_size).possible_neighbors(ri)

        for j in possible_neighbors:
            rj = array([particles[j]['X Pred'],particles[j]['Y Pred'],particles[j]['Z Pred']])
            r_vector = ri - rj
            r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)

            if  r < h:
                # t = 1 -> Fluid - Fluid interaction
                # t = 0 -> Fluid - Boundary interaction
                # t = -1 -> Boundary - Boundary interaction
                if particles[i]['Type'] == 'Water' and particles[j]['Type'] == 'Water':
                    t = 1
                elif (particles[i]['Type'] == 'Water' and particles[j]['Type'] != 'Water'):
                    t = 0
                neighborhood.append([i,j,r_vector,r,t])

        # assigning arrays with interaction type
        fluid_fluid = []
        fluid_boundary = []
        for num, i in enumerate(neighborhood):
            if i[-1] == 1:
                fluid_fluid.append(num)
            elif i[-1] == 0 and neighborhood[i[0]]['Type'] == 'Water':
                fluid_boundary.append(num)

        # Cleaning density and pressure force
        for i in fluid_array:
            particles[i]['Density'] = 0
            particles[neighborhood[i][0]]['Pressure Force']

        ## Calculating predicted density:
        #  For fluid-fluid interactions:
        
        for i in fluid_fluid:
            current_particle = particles[neighborhood[i][0]]
            neighbor = particles[neighborhood[i][1]]
            r_vector = neighborhood[i][2]
            r = neighborhood[i][3]
            t = neighborhood[i][4]
            particles[neighborhood[i][0]]['Density'] += force_fields.Density(current_particle,neighbor,r,h,t,r_vector,"Poly_6")

        #  For fluid-boundary interactions: 
        for i in fluid_boundary:      
            current_particle = particles[neighborhood[i][0]]
            neighbor = particles[neighborhood[i][1]]
            r_vector = neighborhood[i][2]
            r = neighborhood[i][3]
            t = neighborhood[i][4]
            particles[neighborhood[i][0]]['Density'] += force_fields.Density(current_particle,neighbor,r,h,t,r_vector,kernel_name="Poly_6")

        #  Getting max density error
        for i in fluid_array:
            if particles[i]['Density']-rho_0 > max_rho_err:
                max_rho_err = abs(particles[i]['Density']-rho_0)

        # Calculating pressure
        pressure_beta = (delta_t**2 * mass**2 * 2/rho_0**2)**-1
        for i in fluid_array:
            particles[i]['Pressure'] += pressure_delta*pressure_beta*(particles[i]['Density']-rho_0)

        # Calculating pressure force between fluid particles
        for i in fluid_fluid:
            if neighborhood[i][0] != neighborhood[i][1]:
                current_particle = particles[neighborhood[i][0]]
                neighbor = particles[neighborhood[i][1]]
                r_vector = neighborhood[i][2]
                r = neighborhood[i][3]
                t = neighborhood[i][4]
                particles[neighborhood[i][0]]['Pressure Force'] += force_fields.Pressure(current_particle,neighbor,r,h,t,r_vector,"Spiky")

        max_force = 0
        # Calculating new total force and getting max_force
        for i in fluid_array:
            particles[i]['Total Force'] = particles[i]['Pressure Force'] + particles[i]['Viscosity Force'] + particles[i]['Surface Tension Force'] + rho_0 * g
            for j in range(0,3):
                if abs(particles[i]['Total Force'][j]) > max_force:
                    max_force = abs(particles[i]['Total Force'][j])
        k+=1

    # Criterias for delta_t increase
    criteria_1 = 0.19*sqrt(h/max_force) > delta_t
    criteria_2 = max_rho_err < 4.5 * density_fluctuation_permited
    criteria_3 = density_fluctuation_permited < 0.9 * density_fluctuation_permited
    try :
        criteria_4 = 0.39*h/max_vel > delta_t
    except:
        criteria_4 = False

    if criteria_1 and criteria_2 and criteria_3 and criteria_4:
        delta_t += 0.2*delta_t
    else:
    # Criterias for delta_t decrease
        criteria_1 = 0.2*sqrt(h/max_force) < delta_t
        criteria_2 = max_rho_err > 5.5 * density_fluctuation_permited
        criteria_3 = density_fluctuation_permited >= density_fluctuation_permited
        try :
            criteria_4 = 0.4*h/max_vel <= delta_t
        except:
            criteria_4 = False

        if criteria_1 or criteria_2 or criteria_3 or criteria_4:
            delta_t -= 0.2*delta_t
    
    # Calculating all new positions for each particle i
    for i in fluid_array:
        count = 0
        for j in ['X','Y','Z']:
            particles[i][j + " Velocity"] = particles[i][j + " Velocity"] + delta_t * particles[i]['Total Force'][count]/particles[i]['Density']
            particles[i][j] = particles[i][j] + delta_t * particles[i][j + ' Velocity']
            count = count + 1

    time = time + delta_t    

    iteration = iteration + 1
    correction = False
    if iteration % w == 0:
        correction = True
    
    #Saving iterations
    utilities.save_csv(paths[2],iteration,{key: particles[key] for key in fluid_array})
    utilities.save_moving_vtk(paths[0],iteration,{key: particles[key] for key in fluid_array})
    #Making vtk group
    utilities.add_to_group(paths[0],iteration,time,paths[1])
    with open(paths[0].replace("/vtk","/time.txt"), "a") as f:
        f.write(str(time) + "\n")
    utilities.info(time,final_time,start_time,tm.time(),delta_t,iteration)

utilities.save_group(paths[1])