import kernels

# Smoothed density
def Density(current_particle,neighbor_particles,h):
    density = 0
    string_array = ['rx','ry','rz']
    for j in range(0,neighbor_particles.shape[0]):
        for string in string_array:
            #density = density + neighbor_particles.iloc[j]['Mass']*(current_particle['X Velocity']-neighbor_particles.iloc[j]['X Velocity'])*\
            #kernels.Cubic_Spline(abs(neighbor_particles.iloc[j]['rx']),h).Gradient()
            density = density + neighbor_particles.iloc[j]['Mass']*kernels.Cubic_Spline(neighbor_particles.iloc[j][string],h).Kernel()
        
    return density

# Smoothed pressure force
def Pressure(current_particle,neighbor_particles,h,axis):
    
    pressure = 0
    R = 0.006 # Tensile instability term
    string = 'r' + axis.lower()
    for j in range(0,neighbor_particles.shape[0]):
        fac = 1
        if neighbor_particles.iloc[j][string] < 0:
            fac = -1
        pressure = fac * pressure + neighbor_particles.iloc[j]['Mass']* \
        ((current_particle['Pressure']/current_particle['Density']**2)+(neighbor_particles.iloc[j]['Pressure']/neighbor_particles.iloc[j]['Density']**2))* \
        (1+R*(kernels.Cubic_Spline(abs(neighbor_particles.iloc[j][string]),h).Kernel()/kernels.Cubic_Spline(0,h).Kernel())**4)* \
        kernels.Cubic_Spline(neighbor_particles.iloc[j][string],h).Gradient()
    return -pressure

# Smoothed Viscosity

def Viscosity(current_particle,neighbor_particles,h,axis):
    
    viscosity = 0
    string = ['r' + axis.lower(),axis.upper() + ' Velocity']
    for j in range(0,neighbor_particles.shape[0]):
        viscosity = viscosity + neighbor_particles.iloc[j]['Mass']*(neighbor_particles.iloc[j][string[1]]-current_particle[string[1]])/neighbor_particles.iloc[j]['Density']*\
        kernels.Cubic_Spline(neighbor_particles.iloc[j][string[0]],h).Laplacian()
        
    return viscosity

# Smoothed Repulsion Force
def Repulsion(current_particle,neighbor_particles,h,axis):
    repulsion = 0
    string = 'r' + axis.lower()
    for j in range(0,neighbor_particles.shape[0]):
        repulsion = repulsion + neighbor_particles.iloc[j]['Mass']/(neighbor_particles.iloc[j]['Mass']+current_particle['Mass'])*\
        kernels.Cubic_Spline(neighbor_particles.iloc[j][string],h).Kernel()

    return repulsion