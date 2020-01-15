import kernels

# Smoothed density
def Density(current_particle,neighbor_particles,h):
    
    density = 0
    for j in range(0,neighbor_particles.shape[0]):
        #density = density + neighbor_particles.iloc[j]['Mass']*(current_particle['X Velocity']-neighbor_particles.iloc[j]['X Velocity'])*\
        #kernels.Cubic_Spline(abs(neighbor_particles.iloc[j]['rx']),h).Gradient()
        density = density + neighbor_particles.iloc[j]['Mass']*kernels.Cubic_Spline(abs(neighbor_particles.iloc[j]['rx']),h).Gradient()
    
    for j in range(0,neighbor_particles.shape[0]):
        # density = density + neighbor_particles.iloc[j]['Mass']*(current_particle['Y Velocity']-neighbor_particles.iloc[j]['Y Velocity'])*\
        # kernels.Cubic_Spline(abs(neighbor_particles.iloc[j]['ry']),h).Gradient()
        density = density + neighbor_particles.iloc[j]['Mass']*kernels.Cubic_Spline(abs(neighbor_particles.iloc[j]['ry']),h).Gradient()
    
    for j in range(0,neighbor_particles.shape[0]):
        # density = density + neighbor_particles.iloc[j]['Mass']*(current_particle['Z Velocity']-neighbor_particles.iloc[j]['Z Velocity'])*\
        # kernels.Cubic_Spline(abs(neighbor_particles.iloc[j]['rz']),h).Gradient()
        density = density + neighbor_particles.iloc[j]['Mass']*kernels.Cubic_Spline(abs(neighbor_particles.iloc[j]['rz']),h).Gradient()
        
    return density

# Smoothed pressure
def Pressure(current_particle,neighbor_particles,h):
    
    pressure = 0
    for j in range(0,neighbor_particles.shape[0]):
        pressure = pressure + neighbor_particles.iloc[j]['Mass']*(current_particle['Pressure']+neighbor_particles.iloc[j]['Pressure'])/(2*neighbor_particles.iloc[j]['Density'])*\
        kernels.Cubic_Spline(abs(neighbor_particles.iloc[j]['rx']),h).Gradient()
    
    for j in range(0,neighbor_particles.shape[0]):
        pressure = pressure + neighbor_particles.iloc[j]['Mass']*(current_particle['Pressure']+neighbor_particles.iloc[j]['Pressure'])/(2*neighbor_particles.iloc[j]['Density'])*\
        kernels.Cubic_Spline(abs(neighbor_particles.iloc[j]['ry']),h).Gradient()
    
    for j in range(0,neighbor_particles.shape[0]):
        pressure = pressure + neighbor_particles.iloc[j]['Mass']*(current_particle['Pressure']+neighbor_particles.iloc[j]['Pressure'])/(2*neighbor_particles.iloc[j]['Density'])*\
        kernels.Cubic_Spline(abs(neighbor_particles.iloc[j]['rz']),h).Gradient()
    
    return -pressure

def Viscosity(current_particle,neighbor_particles,h):
    
    viscosity = 0
    for j in range(0,neighbor_particles.shape[0]):
        viscosity = viscosity + neighbor_particles.iloc[j]['Mass']*(neighbor_particles.iloc[j]['X Velocity']-current_particle['X Velocity'])/neighbor_particles.iloc[j]['Density']*\
        kernels.Cubic_Spline(abs(neighbor_particles.iloc[j]['rx']),h).Laplacian()
    
    for j in range(0,neighbor_particles.shape[0]):
        viscosity = viscosity + neighbor_particles.iloc[j]['Mass']*(neighbor_particles.iloc[j]['Y Velocity']-current_particle['Y Velocity'])/neighbor_particles.iloc[j]['Density']*\
        kernels.Cubic_Spline(abs(neighbor_particles.iloc[j]['ry']),h).Laplacian()
    
    for j in range(0,neighbor_particles.shape[0]):
        viscosity = viscosity + neighbor_particles.iloc[j]['Mass']*(neighbor_particles.iloc[j]['Z Velocity']-current_particle['Z Velocity'])/neighbor_particles.iloc[j]['Density']*\
        kernels.Cubic_Spline(abs(neighbor_particles.iloc[j]['rz']),h).Laplacian()
        
    return viscosity
            