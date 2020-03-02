
from numpy import sqrt,dot, array, pi, around, arange, mean
import force_fields
import kernels
import utilities
import particle_maker
import time as tm
import hashing

# # Initialization

fluid_radius = 0.01 # Radius of each particle
rho_0 = 1000 # Density at rest 
x = 40 # Number of particles inside the smoothing length
h = round(x**(1/3) * fluid_radius,4)
mass = 4/3*pi*fluid_radius**3*rho_0
# h = 4*fluid_radius # Condition assumed in https://github.com/DeveloperLee/pcisph
boundary_radius = round(h/4,4)
boundary_diameter = boundary_radius*2

# Continuing from last simulation
continue_from_last = False
if continue_from_last == True:
    pull_last = utilities.continue_last_sim()
    particles = pull_last[0]
    particles = particle_maker.make_box([-0.1,-0.1,-0.1],[.41,.41,.41],boundary_radius,0,rho_0,'Stainless Steel',dict_index=len(particles),box=particles)
    iteration = pull_last[1] # Last iteration number
    with open(pull_last[2] + "/time.txt", "r") as f:
        lines = f.read().splitlines()
        time = float(lines[-1])

else:
    # Making particles
    iteration = 1 #First iteration number
    time = 0 # initial time

    # Fluid Particles
    particles = particle_maker.make_prism(0.2,.2,.2,7,7,7,fluid_radius,mass,rho_0,'Water')

    # Boundary Particles
    particles = particle_maker.make_box([0.05,0.05,0.05],[0.5,0.5,0.5],boundary_radius,0,0,'Stainless Steel',dict_index=len(particles),box=particles)

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

#Hashing
hash_table_size = hashing.nextPrime(P*P)
for i in range(0,P):
    point = [particles[i]['X'],particles[i]['Y'],particles[i]['Z']]
    hashing.Hashing(h,hash_table_size)._add(point,i)

# hash_table = hashing.Hashing(h,hash_table_size).d

# Getting neighbors:
boundary_boundary = []
for i in boundary_array:
    ri = array([particles[i]['X'],particles[i]['Y'],particles[i]['Z']])

    possible_neighbors = hashing.Hashing(h,hash_table_size).possible_neighbors(ri)

    for j in possible_neighbors:
        if particles[j]['Type'] != 'Water':
            rj = array([particles[j]['X'],particles[j]['Y'],particles[j]['Z']])
            r_vector = ri - rj
            r = around(sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2),6)

            if  r < h:
                boundary_boundary.append([i,j,r_vector,r,-1])

for i in boundary_array:
    particles[i]['Normal'] = array([0,0,0]).astype(float)

# Calculating the mass of each boundary particle.
for i in boundary_boundary:
    r = around(i[3],10)
    W = kernels.Poly_6(r,h).Kernel()
    particles[i[0]]['Mass'] += W

for i in boundary_array:
    particles[i]['Mass'] = rho_0/particles[i]['Mass']

# Getting boundary max and min positions
min_max_array = array([float('inf'),-float('inf'),float('inf'),-float('inf'),float('inf'),-float('inf')])
for i in boundary_array:
    # max and min array -> [xmin,xmax,ymin,ymax,zmin,zmax]
    x = particles[i]['X']
    y = particles[i]['Y']
    z = particles[i]['Z']

    if x < min_max_array[0]:
        min_max_array[0] = x
    if x > min_max_array[1]:
        min_max_array[1] = x

    if y < min_max_array[2]:
        min_max_array[2] = y
    if y > min_max_array[3]:
        min_max_array[3] = y

    if z < min_max_array[4]:
        min_max_array[4] = z
    if z > min_max_array[5]:
        min_max_array[5] = z

# Calculating the normal vector for each boundary particle
for i in boundary_array:
    normal = array([0.,0.,0.])
    position_vector = [particles[i]['X'],particles[i]['X'],particles[i]['Y'],particles[i]['Y'],particles[i]['Z'],particles[i]['Z']]
    tmp = position_vector/min_max_array
    for j in range(0,len(tmp)):
        if tmp[j] == 1.:
            if j == 0:
                normal[0] = 1.
            elif j == 1:
                normal[0] = -1.
            elif j == 2:
                normal[1] = 1.
            elif j == 3:
                normal[1] = -1.
            elif j == 4:
                normal[2] = 1.
            elif j == 5:
                normal[2] = -1.
    particles[i]['Normal'] = normal


for i in boundary_array:
    particles[i]['Normal'] = particles[i]['Normal']/max(abs(particles[i]['Normal']))

# Used constants
mu = 0.0010518 # Viscosity coefficient
delta_st = 0.0728 # Surface tension coefficient
g = array([0,-9.81,0])
delta_friction = 0.89e-3  # Friction coefficient between 2 interaction surfaces
final_time = 10 #in seconds
episilon = 0.2 #collision dumping function

start_time = tm.time() # actual clock time of initialization

#PCISPH variables
density_fluctuation_permited = rho_0/100 # 1%
maximum_volume_fluctuation = 10*density_fluctuation_permited
delta_shock = 8*density_fluctuation_permited
max_rho_err = float('inf')

max_vel = 0
for i in fluid_array:
    for j in ['X Velocity','Y Velocity','Z Velocity']:
        if abs(particles[i][j]) > max_vel:
            max_vel = abs(particles[i][j])

if max_vel == 0:
    delta_t = 1e-3 # initial delta_t 
else:
    delta_t = 0.25*h/max_vel

# Calculating pressure delta value
project_particle_array = []

for i in arange(-h-fluid_radius,h + 2*fluid_radius,2*fluid_radius):
    for j in arange(-h-fluid_radius,h + 2*fluid_radius,2*fluid_radius):
        for k in arange(-h-fluid_radius,h + 2*fluid_radius,2*fluid_radius):
            project_particle_array.append([i,j,k])

Grad_W = array([0.,0.,0.])
dot_Grad_W = 0
for i in array(project_particle_array):
    r_vector = i
    r = sqrt(i[0]**2 + i[1]**2 + i[2]**2 )
    Grad_W += kernels.Poly_6(r,h).Gradient(r_vector)
    dot_Grad_W += dot(Grad_W,Grad_W) 

pressure_delta = -(dot(Grad_W,Grad_W)-dot_Grad_W)**-1

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
        elif i[-1] == 0 and particles[i[0]]['Type'] == 'Water':
            fluid_boundary.append(num)
    
    # Cleaning Force fields and predicted velocites and positions
    for i in fluid_array:
        particles[i]['Pressure Force'] = array([0.,0.,0.])
        particles[i]['Viscosity Force'] = array([0.,0.,0.])
        particles[i]['Surface Tension Force'] = array([0.,0.,0.])
        particles[i]['Pressure'] = 0
        particles[i]['Normal'] = array([0.,0.,0.])
    
    for i in boundary_array:
        particles[i]['Pressure'] = 0

    # Calculating the surface normal for each particle
    for i in fluid_fluid:
        if neighborhood[i][0] != neighborhood[i][1]:
            current_particle = particles[neighborhood[i][0]]
            neighbor = particles[neighborhood[i][1]]
            r_vector = neighborhood[i][2]
            r = neighborhood[i][3]
            Grad_W = kernels.Poly_6(r,h).Gradient(r_vector)
            particles[neighborhood[i][0]]['Normal'] += around(h*neighbor['Mass']/neighbor['Density'] * Grad_W,5)

    for i in fluid_boundary:
        if neighborhood[i][0] != neighborhood[i][1]:
            current_particle = particles[neighborhood[i][0]]
            neighbor = particles[neighborhood[i][1]]
            r_vector = neighborhood[i][2]
            r = neighborhood[i][3]
            Grad_W = kernels.Poly_6(r,h).Gradient(r_vector)
            particles[neighborhood[i][0]]['Normal'] += around(h*neighbor['Mass']/neighbor['Density'] * Grad_W)

    # Computing viscosity and surface tension forces    
    for i in fluid_fluid:
        if neighborhood[i][0] != neighborhood[i][1]:
            current_particle = particles[neighborhood[i][0]]
            neighbor = particles[neighborhood[i][1]]
            r_vector = neighborhood[i][2]
            r = neighborhood[i][3]
            t = neighborhood[i][4]
            particles[neighborhood[i][0]]['Viscosity Force'] += force_fields.Viscosity(current_particle,neighbor,r,h,t,"Viscosity_Kernel",mu)
            particles[neighborhood[i][0]]['Surface Tension Force'] += force_fields.Surface_Tension(current_particle,neighbor,r,h,r_vector,"Poly_6",delta_st,rho_0)

    for i in fluid_boundary:
        if neighborhood[i][0] != neighborhood[i][1]:
            current_particle = particles[neighborhood[i][0]]
            neighbor = particles[neighborhood[i][1]]
            r_vector = neighborhood[i][2]
            r = neighborhood[i][3]
            t = neighborhood[i][4]
            particles[neighborhood[i][0]]['Viscosity Force'] += force_fields.Viscosity(current_particle,neighbor,r,h,t,"Viscosity_Kernel",mu)
            particles[neighborhood[i][0]]['Surface Tension Force'] += force_fields.Surface_Tension(current_particle,neighbor,r,h,r_vector,"Poly_6",delta_st,rho_0)
    
    # Calculating summation of all non-pressure forces
    for i in fluid_array:
        particles[i]['Total Force'] = particles[i]['Viscosity Force'] + particles[i]['Surface Tension Force'] + g*particles[i]['Mass']

    max_rho_err_t_minus_1 = max_rho_err

    k=0
    while k < 3:
        
        # Predicting velocity and position (t+1)
        for i in fluid_array:
            count = 0
            for j in ['X','Y','Z']:
                particles[i][j + " Pred Velocity"] = particles[i][j + " Velocity"] + delta_t * particles[i]['Total Force'][count]/particles[i]['Mass']
                particles[i][j + " Pred"] = particles[i][j] + delta_t * particles[i][j + ' Pred Velocity']
                count = count + 1

        # New neighborhood calculation with predicted positions
        # Hashing
        for i in fluid_array:
            point = [particles[i]['X Pred'],particles[i]['Y Pred'],particles[i]['Z Pred']]
            hashing.Hashing(h,hash_table_size)._add(point,i)

        # neighborhood search
        neighborhood = []
        colliders = []
        for i in fluid_array:
            ri = array([particles[i]['X Pred'],particles[i]['Y Pred'],particles[i]['Z Pred']])

            possible_neighbors = hashing.Hashing(h,hash_table_size).possible_neighbors(ri)

            for j in possible_neighbors:
                if particles[j]['Type'] != 'Water':
                    rj = array([particles[j]['X'],particles[j]['Y'],particles[j]['Z']])
                else:
                    rj = array([particles[j]['X Pred'],particles[j]['Y Pred'],particles[j]['Z Pred']])
                r_vector = ri - rj
                r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)

                if  r < h:
                    # t = 1 -> Fluid - Fluid interaction
                    # t = 0 -> Fluid - Boundary interaction
                    if particles[i]['Type'] == 'Water' and particles[j]['Type'] == 'Water':
                        t = 1
                    elif particles[i]['Type'] == 'Water' and particles[j]['Type'] != 'Water':
                        t = 0
                    neighborhood.append([i,j,r_vector,r,t])
                
                # Predicting world collision
                if r < boundary_diameter:
                    if particles[i]['Type'] == 'Water' and particles[j]['Type'] != 'Water':
                        #colliders -> [fluid particle , boundary particle, r_vector, r]
                        colliders.append([i,j,r_vector,r])

        # Assigning arrays with interaction type
        fluid_fluid = []
        fluid_boundary = []
        for num, i in enumerate(neighborhood):
            if i[-1] == 1:
                fluid_fluid.append(num)
            elif i[-1] == 0 and particles[i[0]]['Type'] == 'Water':
                fluid_boundary.append(num)

        # Cleaning density and pressure force
        for i in fluid_array:
            particles[i]['Density'] = 0
            particles[i]['Pressure Force'] = array([0.,0.,0.])

        for i in boundary_array:
            particles[i]['Density'] = 0
            particles[i]['Pressure Force'] = array([0.,0.,0.])

        ## Calculating predicted density:
        #  For fluid-fluid interactions:
        
        for i in fluid_fluid:
            current_particle = particles[neighborhood[i][0]]
            neighbor = particles[neighborhood[i][1]]
            r_vector = neighborhood[i][2]
            r = neighborhood[i][3]
            particles[neighborhood[i][0]]['Density'] += force_fields.Density(current_particle,neighbor,r,h,r_vector,"Poly_6")

        #  For fluid-boundary interactions: 
        for i in fluid_boundary:      
            current_particle = particles[neighborhood[i][0]]
            neighbor = particles[neighborhood[i][1]]
            r_vector = neighborhood[i][2]
            r = neighborhood[i][3]
            density = force_fields.Density(current_particle,neighbor,r,h,r_vector,"Poly_6")
            particles[neighborhood[i][0]]['Density'] += density
            particles[neighborhood[i][1]]['Density'] += density

        for i in boundary_array:
            particles[i]['Density'] += 1000

        #  Getting max density error
        max_rho_err = 0
        density_error_array = []
        for i in fluid_array:
            density_error_array.append(max(particles[i]['Density']-rho_0,0))
            if particles[i]['Density']-rho_0 > max_rho_err:
                max_rho_err = max(particles[i]['Density']-rho_0,0)

        for i in boundary_array:
            density_error_array.append(max(particles[i]['Density']-rho_0,0))
            if particles[i]['Density']-rho_0 > max_rho_err:
                max_rho_err = max(particles[i]['Density']-rho_0,0)

        # Calculating pressure
        pressure_beta = (delta_t**2 * mass**2 * 2/rho_0**2)**-1
        for i in fluid_array:
            particles[i]['Pressure'] += max(pressure_delta*pressure_beta*(particles[i]['Density']-rho_0),0)
 
        for i in boundary_array:
            particles[i]['Pressure'] += max(pressure_delta*pressure_beta*(particles[i]['Density']-rho_0),0)

        # if max_rho_err > 0:
        #     print()
        # Calculating pressure force between fluid particles
        for i in fluid_fluid:
            if neighborhood[i][0] != neighborhood[i][1]:
                current_particle = particles[neighborhood[i][0]]
                neighbor = particles[neighborhood[i][1]]
                r_vector = neighborhood[i][2]
                r = neighborhood[i][3]
                t = neighborhood[i][4]
                particles[neighborhood[i][0]]['Pressure Force'] += force_fields.Pressure(current_particle,neighbor,r,h,t,r_vector,"Spiky")

        # Calculating pressure forces between fluid and particle interactions
        for i in fluid_boundary:
            current_particle = particles[neighborhood[i][0]]
            neighbor = particles[neighborhood[i][1]]
            r_vector = neighborhood[i][2]
            r = neighborhood[i][3]
            t = neighborhood[i][4]
            particles[neighborhood[i][0]]['Pressure Force'] += force_fields.Pressure(current_particle,neighbor,r,h,t,r_vector,"Spiky")

        max_force = 0
        # Calculating new total force and getting max_force
        for i in fluid_array:
            particles[i]['Total Force'] = particles[i]['Pressure Force'] + particles[i]['Viscosity Force'] + particles[i]['Surface Tension Force'] + g*particles[i]['Mass']
            for j in range(0,3):
                if abs(particles[i]['Total Force'][j]) > max_force:
                    max_force = abs(particles[i]['Total Force'][j])
        k+=1
        # print(max_rho_err)
    
    # Calculating all new velocities and positions for each particle i
    for i in fluid_array:
        count = 0
        for j in ['X','Y','Z']:
            particles[i][j + " Velocity"] = particles[i][j + " Velocity"] + delta_t * particles[i]['Total Force'][count]/particles[i]['Mass']
            particles[i][j] = particles[i][j] + delta_t * particles[i][j + ' Velocity']
            count = count + 1
    
    # Computing world collision
    count = 0
    for a in colliders:
        i = a[0]
        j = a[1]
        r = a[3]
        wib = (boundary_diameter-r)/boundary_diameter
        wib2 = wib*(boundary_diameter-r)
        ni = wib*particles[j]['Normal']
        colliders[count].append(ni)
        colliders[count].append(wib)
        colliders[count].append(wib2)
        count += 1

    for a in colliders:
        i = a[0]
        particles[i]['ni'] = array([0.,0.,0.])
        particles[i]['wib'] = 0
        particles[i]['wib2'] = 0
    
    for a in colliders:
        i = a[0]
        particles[i]['ni'] += a[4]
        particles[i]['wib'] += a[5]
        particles[i]['wib2'] += a[6]
    
    computed = []
    for a in colliders:
        i = a[0]
        if i not in computed:
            j = a[1]
            r = a[3]
            wib = (boundary_diameter-r)/boundary_diameter
            AVG_boundary_normals_abs = sqrt(dot(particles[i]['ni'],particles[i]['ni']))
            v = [particles[i]["X Velocity"],particles[i]["Y Velocity"],particles[i]["Z Velocity"]]
            n = particles[j]['Normal']
            vt = dot(dot(v,n),n)
            count = 0
            for axis in ['X','Y','Z']:
                particles[i][axis] += particles[i]['wib2']/particles[i]['wib']*(particles[i]['ni'][count]/AVG_boundary_normals_abs)
                particles[i][axis + " Velocity"] = episilon*(particles[i][axis + " Velocity"] - vt[count])
                count = count + 1
            computed.append(i)
        try:
            particles[i].pop('ni')
            particles[i].pop('wib')
            particles[i].pop('wib2')
        except:
            pass

    # Getting max velocity to calculate delta t
    max_vel = 0
    for i in fluid_array:
        for j in ['X Velocity','Y Velocity','Z Velocity']:
            if abs(particles[i][j]) > max_vel:
                max_vel = abs(particles[i][j])

    time = time + delta_t    

    iteration = iteration + 1

    # Shock detection
    criteria_1 = max_rho_err - max_rho_err_t_minus_1 > delta_shock
    criteria_2 = max_rho_err > maximum_volume_fluctuation
    criteria_3 = 0.45*h/max_vel < delta_t

    if criteria_1 or criteria_2 or criteria_3:
        print("Shock detected! Going 2 iterations backwards!")
        delta_t = min(0.2*sqrt(h/max_force),0.25*h/max_vel)
        iteration -= 2
        

    # Criterias for delta_t increase
    try:
        criteria_1 = 0.19*sqrt(h/max_force) > delta_t
    except:
        criteria_1 = False
    criteria_2 = max_rho_err < 4.5 * density_fluctuation_permited
    criteria_3 = mean(density_error_array) < 0.9 * density_fluctuation_permited
    try :
        criteria_4 = 0.39*h/max_vel > delta_t
    except:
        criteria_4 = False

    if criteria_1 and criteria_2 and criteria_3 and criteria_4:
        delta_t += 0.002*delta_t
    else:
    # Criterias for delta_t decrease
        try:
            criteria_1 = 0.2*sqrt(h/max_force) < delta_t
        except:
            criteria_1 = False
        criteria_2 = max_rho_err > 5.5 * density_fluctuation_permited
        criteria_3 = mean(density_error_array) >= density_fluctuation_permited
        try :
            criteria_4 = 0.4*h/max_vel <= delta_t
        except:
            criteria_4 = False

        if criteria_1 or criteria_2 or criteria_3 or criteria_4:
            delta_t -= 0.002*delta_t
    
    if time + delta_t > final_time:
        delta_t = final_time - time

    #Saving iterations
    utilities.save_csv(paths[2],iteration,{key: particles[key] for key in fluid_array})
    utilities.save_moving_vtk(paths[0],iteration,{key: particles[key] for key in fluid_array})
    #Making vtk group
    utilities.add_to_group(paths[0],iteration,time,paths[1])
    with open(paths[0].replace("/vtk","/time.txt"), "a") as f:
        f.write(str(time) + "\n")
    utilities.info(time,final_time,start_time,tm.time(),delta_t,iteration,max_rho_err)

utilities.save_group(paths[1])