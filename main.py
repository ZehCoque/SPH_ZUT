
from numpy import sqrt, array, pi, dot, sign, around
import csv
import force_fields
import kernels
import utilities
import particle_maker
import time as tm

# # Initialization

# Making particles

fluid_radius = 0.01 # Radius of each particle
boundary_radius = 0.01
rho_0 = 1000 # Density at rest
mass = 1 # temporary mass of each particle

boundary = particle_maker.make_box([-0.1,-0.1,-0.1],[.3,.3,.3],boundary_radius,0,'Stainless Steel',dict_index=0,box={})

moving = particle_maker.make_prism2([0,0,0],[0.1,.2,.1],fluid_radius,mass,'Water',dict_index=0,prism={})

N = len(moving)

print()
print("#"*40)
print('Initializing with:')
print(str(N) + ' fluid particles and ' + str(len(boundary)) + ' boundary particles')

#Calculating Kernel Support radius
V = N * fluid_radius**3 * 4/3 * pi
h = (3*V*25/(4*pi*N))**(1/3)

# Calculating desired mass-density of each particle:
print("#"*40)
print('Calculating desired mass-density of each particle')
# We want the fluid in the initial conguration to exist roughly at the reference density. 
# One way to do this is to take the volume in the indicated body of fuid, multiply by the mass density, and divide by the
# number of particles; but that requires that we be able to compute the volume of the fluid region. 
# Alternately, we can simply compute the average mass density assuming each particle has mass one, then use that to compute the particle
# mass necessary in order to achieve the desired reference density. Applications of Parallel Computers (Bindel, 2011)

rho = []
rhos = 0
rhos2 = 0

for i in range(0,N):
    moving[i]['Fluid Neighbors'] = []
    moving[i]['moving_r'] = []
    ri_x = moving[i]['X']
    ri_y = moving[i]['Y']
    ri_z = moving[i]['Z']

    for j in range(0,N):
            rx = ri_x-moving[j]['X']
            ry = ri_y-moving[j]['Y']
            rz = ri_z-moving[j]['Z']

            if sqrt(rx**2+ry**2+rz**2) <= h and i!=j:
                moving[i]['moving_r'].append([rx,ry,rz])
                moving[i]['Fluid Neighbors'].append(j)

for i in range(0,N):
    rho.append(force_fields.Density(moving,i,h,True,rho_0,'Poly_6'))
    rhos += rho[i]
    rhos2 += rho[i]**2

mass = rho_0*rhos/rhos2
for i in range(0,N):
    moving[i]['Mass'] = mass

# Calculating the volume of each boundary particle
print("#"*40)
print('Calculating the volume of each boundary particle')
for i in range(0,len(boundary)):
    ri_x = boundary[i]['X']
    ri_y = boundary[i]['Y']
    ri_z = boundary[i]['Z']
    W = 0
    boundary[i]['psi'] = 0
    for j in range(0,len(boundary)):
        rx = ri_x-boundary[j]['X']
        ry = ri_y-boundary[j]['Y']
        rz = ri_z-boundary[j]['Z']
        r = sqrt(rx**2+ry**2+rz**2)

        if r <= h:
            W += kernels.Poly_6(r,h).Kernel()
    boundary[i]['psi'] = rho_0/W

#  -Gets initial boundary neighbors for each moving particle i
print("#"*40)
print('Getting initial boundary neighbors for each moving particle i')
for i in range(0,N):
    moving[i]['Boundary Neighbors'] = []
    moving[i]['boundary_r'] = []
    ri_x = moving[i]['X']
    ri_y = moving[i]['Y']
    ri_z = moving[i]['Z']
    for j in range(0,len(boundary)):
            rx = ri_x-boundary[j]['X']
            ry = ri_y-boundary[j]['Y']
            rz = ri_z-boundary[j]['Z']

            if sqrt(rx**2+ry**2+rz**2) <= h:
                moving[i]['boundary_r'].append([rx,ry,rz])
                moving[i]['Boundary Neighbors'].append(j)

# Used constants
k = 3 # Gas stiffness
mu = 3.5*1.003e-3 # Viscosity coefficient
delta = 0.0728 # Surface tension coefficient
g = array([0,-9.81,0])
lam_v = 0.4 # delta_t for velocity parameter
lam_f = 0.2 # delta_t for force parameter
gamma = 7 # power of pressure calculation
time = 0 # initial time
final_time = 1 #in seconds
iteration = 1 #First iteration number
start_time = tm.time() # actual clock time of initialization
alpha = 0.01 # First artificial viscosity term
beta = 0 # Second artificial viscosity term
w = 50 # steps for density correction
correction = True # True for correction in the first iteration
delta_friction = 1 # Friction coefficient between 2 interaction surfaces
c_R = 0.01 # coefficient of restitution (between 0 and 1)

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
utilities.save_csv(paths[2],iteration,moving)
utilities.save_moving_vtk(paths[0],iteration,moving)
utilities.save_boundary_vtk(paths[0],boundary)
utilities.add_to_group(paths[0],iteration,time,paths[1])

print("#"*40)

# Stop when simulation time reaches final time
while time < final_time:
    
    # First iteration through all particles
    #  -Gets moving neighbors for each moving particle i
    for i in range(0,N):
        moving[i]['Fluid Neighbors'] = []
        moving[i]['moving_r'] = []
        ri_x = moving[i]['X']
        ri_y = moving[i]['Y']
        ri_z = moving[i]['Z']

        for j in range(0,N):
                rx = ri_x-moving[j]['X']
                ry = ri_y-moving[j]['Y']
                rz = ri_z-moving[j]['Z']

                if sqrt(rx**2+ry**2+rz**2) <= h and i!=j:
                    moving[i]['moving_r'].append([rx,ry,rz])
                    moving[i]['Fluid Neighbors'].append(j)

    # Second iteration through all particles
    #  -Calculates density and pressure for each particle i
    for i in range(0,N):
        moving[i]['Density'] = force_fields.Density(moving,i,h,correction,rho_0,'Poly_6',boundary)
        moving[i]['Pressure'] = k*((moving[i]['Density']/rho_0)**gamma-1)

    # Third iteration through all particles
    #  -Calculates all INTERNAL force fields (pressure,viscosity and others) for each particle i
    for i in range(0,N):
        moving[i]['Pressure Force'] = force_fields.Pressure(moving,i,h,alpha,beta,c,'Spiky')
        moving[i]['Viscosity Force'] = force_fields.Viscosity_Kernel(moving,i,h,mu,'Viscosity')
        

    # Fourth iteration through all particles
    #  -Calculates all EXTERNAL force fields (pressure,viscosity and others) for each particle i
    for i in range(0,N):
        moving[i]['Surface Tension Force'] = force_fields.Surface_Tension(moving,i,h,delta,'Poly_6')
        neighbors = { key: boundary[key] for key in moving[i]['Boundary Neighbors'] }
        moving[i]['Boundary-Fluid Pressure'] = 0
        moving[i]['Boundary-Fluid Friction'] = 0
        if len(neighbors) > 0:
            count = 0
            for j in neighbors:
                moving[i]['Boundary-Fluid Pressure'] += force_fields.Boundary_Fluid_Pressure(moving[i],neighbors[j],h,count,'Poly_6')
                moving[i]['Boundary-Fluid Friction'] += force_fields.Boundary_Fluid_Friction(moving[i],neighbors[j],h,delta_friction,c,count,'Poly_6')
                count += 1

    max_force = 0
    # Calculating total force for each axis
    for i in range(0,N):
            moving[i]['Total Force'] = moving[i]['Pressure Force'] + moving[i]['Viscosity Force'] + moving[i]['Surface Tension Force'] + rho_0 * g + moving[i]['Boundary-Fluid Pressure'] + moving[i]['Boundary-Fluid Friction']
            for j in range(0,3):
                if abs(moving[i]['Total Force'][j]) > abs(max_force):
                    max_force = abs(moving[i]['Total Force'][j])
                    
    max_vel = 0
    for i in range(0,N):
        for j in ['X Velocity','Y Velocity','Z Velocity']:
            if abs(moving[i][j]) > max_vel:
                max_vel = abs(moving[i][j])

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
    for i in range(0,N):
        count = 0
        for j in ['X','Y','Z']:
            moving[i][j + " Velocity"] = moving[i][j + " Velocity"] + delta_t * moving[i]['Total Force'][count]/moving[i]['Density']
            moving[i][j] = moving[i][j] + delta_t * moving[i][j + ' Velocity']
            count = count + 1

    # Collision Handling

    # Collision detection
    colliders = []
    for i in range(0,N):
        moving[i]['Boundary Neighbors'] = []
        moving[i]['boundary_r'] = []
        ri_x = moving[i]['X']
        ri_y = moving[i]['Y']
        ri_z = moving[i]['Z']
        for j in range(0,len(boundary)):
                rx = ri_x-boundary[j]['X']
                ry = ri_y-boundary[j]['Y']
                rz = ri_z-boundary[j]['Z']
                r = sqrt(rx**2+ry**2+rz**2)
                if r <= h:
                    moving[i]['boundary_r'].append([rx,ry,rz])
                    moving[i]['Boundary Neighbors'].append(j)
                
                if r < (fluid_radius + boundary_radius):
                    colliders.append([i,j])

    # Collision Response
    try:
        current_particle = colliders[0][0]
    except:
        pass
  
    for i in range(0,len(colliders)):

        if current_particle == colliders[i][0]:
            fluid_collider = moving[colliders[i][0]]
            boundary_collider = boundary[colliders[i][1]]

            r_index = fluid_collider['Boundary Neighbors'].index(colliders[i][1])
            string = ['X Velocity','Y Velocity','Z Velocity']
            d = abs(sqrt(fluid_collider['boundary_r'][r_index][0]**2+fluid_collider['boundary_r'][r_index][1]**2+fluid_collider['boundary_r'][r_index][2]**2) - (fluid_radius + boundary_radius))
            current_particle = colliders[i][0]
            new_vel = array([0.,0.,0.])
            mag_old_vel = sqrt(fluid_collider[string[0]]**2+fluid_collider[string[1]]**2+fluid_collider[string[2]]**2)
            direction = sign(around(fluid_collider['boundary_r'][r_index], decimals = 10))
            normal = -1*direction
            
            old_vel = []
            for v_axis in string:
                old_vel.append(fluid_collider[v_axis])
            
            new_vel = old_vel - (1+c_R*d/(delta_t*mag_old_vel))*dot(old_vel,normal)*normal

            count = 0
            for v_axis in string:
                moving[colliders[i][0]][v_axis] = new_vel[count] * lam_v * fluid_radius*2 /max(new_vel)
                moving[colliders[i][0]][v_axis[0]] = moving[colliders[i][0]][v_axis[0]] + lam_v * fluid_radius*2 /max(new_vel) * new_vel[count]
                count += 1
             

    time = time + delta_t    

    iteration = iteration + 1
    correction = False
    if iteration % w == 0:
        correction = True
    
    #Saving each iteration
    utilities.save_csv(paths[2],iteration,moving)
    utilities.save_moving_vtk(paths[0],iteration,moving)
    
    #Making vtk group
    utilities.add_to_group(paths[0],iteration,time,paths[1])

    utilities.info(time,final_time,start_time,tm.time(),delta_t)

utilities.save_group(paths[1])