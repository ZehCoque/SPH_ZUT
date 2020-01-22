import kernels


# Smoothed density
def Density(dict_moving,dict_boundary,i,h):
    density = 0
    neighbors = { key: dict_moving[key] for key in dict_moving[i]['Fluid Neighbors'] }
    count = 0
    for j in dict_moving[i]['Fluid Neighbors']:
        for axis in range(0,3):
            density = density + neighbors[j]['Mass']*kernels.Cubic_Spline(dict_moving[i]['moving_r'][count][axis],h).Kernel()
        count = count + 1
    
    neighbors = { key: dict_boundary[key] for key in dict_moving[i]['Boundary Neighbors'] }
    count = 0
    for j in dict_moving[i]['Boundary Neighbors']:
        for axis in range(0,3):
            density = density + neighbors[j]['Mass']*kernels.Cubic_Spline(dict_moving[i]['boundary_r'][count][axis],h).Kernel()
        count = count + 1
    
    return density

# Smoothed pressure force
def Pressure(dict_moving,dict_boundary,i,h):
    R = 0.006 # Tensile instability term

    pressure = [0.,0.,0.]
    neighbors = { key: dict_moving[key] for key in dict_moving[i]['Fluid Neighbors'] }
    count = 0
    for j in dict_moving[i]['Fluid Neighbors']:
        for axis in range(0,3):
            fac = 1
            if dict_moving[i]['moving_r'][count][axis] > 0:
                fac = -1
            elif dict_moving[i]['moving_r'][count][axis] == 0:
                continue
            pressure[axis] = pressure[axis] - fac * neighbors[j]['Mass']* \
            ((dict_moving[i]['Pressure']/dict_moving[i]['Density']**2)+(neighbors[j]['Pressure']/neighbors[j]['Density']**2))* \
            (1+R*(kernels.Cubic_Spline(dict_moving[i]['moving_r'][count][axis],h).Kernel()/kernels.Cubic_Spline(0,h).Kernel())**4)* \
            kernels.Cubic_Spline(dict_moving[i]['moving_r'][count][axis],h).Gradient()
        count = count + 1

    neighbors = { key: dict_boundary[key] for key in dict_moving[i]['Boundary Neighbors'] }
    count = 0
    for j in dict_moving[i]['Boundary Neighbors']:
        for axis in range(0,3):
            fac = 1
            if dict_moving[i]['boundary_r'][count][axis] > 0:
                fac = -1
            elif dict_moving[i]['boundary_r'][count][axis] == 0:
                continue
            pressure[axis] = pressure[axis] - fac * neighbors[j]['Mass']* \
            2*dict_moving[i]['Pressure']/dict_moving[i]['Density']**2* \
            (1+R*(kernels.Cubic_Spline(dict_moving[i]['boundary_r'][count][axis],h).Kernel()/kernels.Cubic_Spline(0,h).Kernel())**4)* \
            kernels.Cubic_Spline(dict_moving[i]['boundary_r'][count][axis],h).Gradient()
        count = count + 1

    return pressure

# Smoothed Viscosity

def Viscosity(dict_moving,i,h,mu):
    viscosity = [0.,0.,0.]
    neighbors = { key: dict_moving[key] for key in dict_moving[i]['Fluid Neighbors'] }
    count = 0
    for j in dict_moving[i]['Fluid Neighbors']:
        axis_num = 0
        for axis in ['X Velocity','Y Velocity','Z Velocity']:
            viscosity = viscosity + neighbors[j]['Mass']*(neighbors[j][axis]-dict_moving[i][axis])/neighbors[j]['Density']*\
            kernels.Cubic_Spline(dict_moving[i]['moving_r'][count][axis_num],h).Laplacian()
            axis_num = axis_num+1
        count = count + 1
        
    return mu*viscosity

# Smoothed Repulsion Force
def Repulsion(current_particle,neighbor_particles,h,axis):
    repulsion = 0
    string = 'r' + axis.lower()
    for j in range(0,neighbor_particles.shape[0]):
        repulsion = repulsion + neighbor_particles.iloc[j]['Mass']/(neighbor_particles.iloc[j]['Mass']+current_particle['Mass'])*\
        kernels.Cubic_Spline(neighbor_particles.iloc[j][string],h).Kernel()

    return repulsion
