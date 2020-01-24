from kernels import *
from numpy import array, sqrt

# Smoothed density
def Density(dict_moving,dict_boundary,i,h,kernel_name):
    # Initialization
    density = 0
    kernel_name = globals()['%s' % kernel_name]
    neighbors = { key: dict_moving[key] for key in dict_moving[i]['Fluid Neighbors'] }
    count = 0
    
    # Fluid Neighbors
    for j in dict_moving[i]['Fluid Neighbors']:
        r_vector = array(dict_moving[i]['moving_r'][count])
        r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)
        W = kernel_name(r,h).Kernel()

        density = density + neighbors[j]['Mass']*W
    
    # Boundary Neighbors
    neighbors = { key: dict_boundary[key] for key in dict_moving[i]['Boundary Neighbors'] }
    for j in dict_moving[i]['Boundary Neighbors']:
        r_vector = array(dict_moving[i]['boundary_r'][count])
        r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)
        W = kernel_name(r,h).Kernel()

        density = density + neighbors[j]['Mass']*W

    return density + dict_moving[i]['Mass']*kernel_name(0,h).Kernel()

# Smoothed pressure force
def Pressure(dict_moving,dict_boundary,i,h,kernel_name):
    # Initialization
    #R = 0.006 # Tensile instability term
    kernel_name = globals()['%s' % kernel_name]
    pressure = array([0.,0.,0.])
    direction = array([0.,0.,0.])
    count = 0

    # Fluid neighbors
    neighbors = { key: dict_moving[key] for key in dict_moving[i]['Fluid Neighbors'] }
    for j in dict_moving[i]['Fluid Neighbors']:
        r_vector = array(dict_moving[i]['moving_r'][count])
        r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)
        Grad_W = kernel_name(r,h).Gradient()
        for dir_num in range(0,len(direction)):
            if r_vector[dir_num] == 0:
                direction[dir_num] = 0
            else: 
                direction[dir_num] = (r_vector[dir_num])/abs(r_vector[dir_num])

        pressure = pressure + direction * neighbors[j]['Mass']* \
        ((dict_moving[i]['Pressure']/dict_moving[i]['Density']**2)+(neighbors[j]['Pressure']/neighbors[j]['Density']**2))*Grad_W
        count += 1

    # Boundary Neighbors
    neighbors = { key: dict_boundary[key] for key in dict_moving[i]['Boundary Neighbors'] }
    count = 0
    for j in dict_moving[i]['Boundary Neighbors']:
        r_vector = array(dict_moving[i]['boundary_r'][count])
        r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)
        W = kernel_name(r,h).Gradient()
        for dir_num in range(0,len(direction)):
            if r_vector[dir_num] == 0:
                direction[dir_num] = 0
            else: 
                direction[dir_num] = (r_vector[dir_num])/abs(r_vector[dir_num])

        pressure = pressure + direction * neighbors[j]['Mass']*(dict_moving[i]['Pressure']/dict_moving[i]['Density']**2)*W
        count += 1

    return dict_moving[i]['Density']*pressure

# Smoothed Viscosity

def Viscosity_Kernel(dict_moving,i,h,mu,kernel_name):
    #Initialization
    kernel_name = globals()['%s' % kernel_name]
    viscosity = array([0.,0.,0.])

    # Fluid Neighbors
    neighbors = { key: dict_moving[key] for key in dict_moving[i]['Fluid Neighbors'] }
    count = 0
    for j in dict_moving[i]['Fluid Neighbors']:
        r_vector = array(dict_moving[i]['moving_r'][count])
        r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)
        Laplacian_W = kernel_name(r,h).Laplacian()

        dv = array([neighbors[j]['X Velocity']-dict_moving[i]['X Velocity'],neighbors[j]['Y Velocity']-dict_moving[i]['Y Velocity'],neighbors[j]['Z Velocity']-dict_moving[i]['Z Velocity']])

        viscosity = viscosity + neighbors[j]['Mass']*dv/neighbors[j]['Density']*Laplacian_W
        
    return mu*viscosity

# Smoothed Surface Tension Force
def Surface_Tension(dict_moving,i,h,delta,kernel_name):
    # Initialization
    kernel_name = globals()['%s' % kernel_name]
    surface = array([0.,0.,0.])
    direction = array([0.,0.,0.])
    count = 0

    neighbors = { key: dict_moving[key] for key in dict_moving[i]['Fluid Neighbors'] }
    for j in dict_moving[i]['Fluid Neighbors']:
        r_vector = array(dict_moving[i]['moving_r'][count])
        r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)
        Grad_W = kernel_name(r,h).Gradient()
        Laplacian_W = kernel_name(r,h).Laplacian()

        n = neighbors[j]['Mass']/neighbors[j]['Density'] * Grad_W
        Laplacian_c = neighbors[j]['Mass']/neighbors[j]['Density'] * Laplacian_W

        for dir_num in range(0,len(direction)):
            if r_vector[dir_num] == 0:
                direction[dir_num] = 0
            else: 
                direction[dir_num] = (r_vector[dir_num])/abs(r_vector[dir_num])
        
        surface = surface - direction * delta * Laplacian_c * n/abs(n)

    return surface
