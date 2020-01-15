import kernels

# Smoothed density
def Density(current_particle,neighbor_particles,h):
    
    density = 0
    string_array = ['rx','ry','rz']
    for j in range(0,neighbor_particles.shape[0]):
        for string in string_array:
            #density = density + neighbor_particles.iloc[j]['Mass']*(current_particle['X Velocity']-neighbor_particles.iloc[j]['X Velocity'])*\
            #kernels.Cubic_Spline(abs(neighbor_particles.iloc[j]['rx']),h).Gradient()
            density = density + neighbor_particles.iloc[j]['Mass']*kernels.Cubic_Spline(abs(neighbor_particles.iloc[j][string]),h).Gradient()
        
    return density

# Smoothed pressure
def Pressure(current_particle,neighbor_particles,h,axis):
    
    pressure = 0
    string = 'r' + axis.lower()
    for j in range(0,neighbor_particles.shape[0]):
        pressure = pressure + neighbor_particles.iloc[j]['Mass']*(current_particle['Pressure']+neighbor_particles.iloc[j]['Pressure'])/(2*neighbor_particles.iloc[j]['Density'])*\
        kernels.Cubic_Spline(abs(neighbor_particles.iloc[j][string]),h).Gradient()
    
    return -pressure

def Viscosity(current_particle,neighbor_particles,h,axis):
    
    viscosity = 0
    string = ['r' + axis.lower(),axis.upper() + ' Velocity']
    for j in range(0,neighbor_particles.shape[0]):
        viscosity = viscosity + neighbor_particles.iloc[j]['Mass']*(neighbor_particles.iloc[j][string[1]]-current_particle[string[1]])/neighbor_particles.iloc[j]['Density']*\
        kernels.Cubic_Spline(abs(neighbor_particles.iloc[j][string[0]]),h).Laplacian()
        
    return viscosity
            