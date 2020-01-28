from kernels import *
from numpy import array, sqrt, dot, sign, around

# Smoothed density
def Density(dict_moving,i,h,correction,rho_0,kernel_name,dict_boundary={}):
    # Initialization
    density = 0
    kernel_name = globals()['%s' % kernel_name]
    try:
        neighbors = { key: dict_moving[key] for key in dict_moving[i]['Fluid Neighbors'] }
    except:
        dict_moving.pop(i)
        neighbors = dict_moving

    count = 0 
    if correction:
        beta = Kernel_Correction(neighbors,kernel_name,h)
    else:
        beta = []

    # Fluid Neighbors
    for j in neighbors:
        r_vector = array(dict_moving[i]['moving_r'][count])
        r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)
        if correction:
            W = (beta[0]+beta[1]*r_vector[0]+beta[2]*r_vector[1]+beta[3]*r_vector[2])*kernel_name(r,h).Kernel()
        else:
            W = kernel_name(r,h).Kernel()

        density += neighbors[j]['Mass']*W
        count += 1

    count = 0 
    #Boundary Neighbors
    if dict_boundary != {}:
        neighbors = { key: dict_boundary[key] for key in dict_moving[i]['Boundary Neighbors'] }
        for j in neighbors:
            r_vector = array(dict_moving[i]['boundary_r'][count])
            r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)
            if correction:
                W = (beta[0]+beta[1]*r_vector[0]+beta[2]*r_vector[1]+beta[3]*r_vector[2])*kernel_name(r,h).Kernel()
            else:
                W = kernel_name(r,h).Kernel()
            
            density += dict_boundary[j]['psi'] * W
            count += 1

    return density + dict_moving[i]['Mass']*kernel_name(0,h).Kernel()

# Smoothed pressure force
def Pressure(dict_moving,i,h,alpha,beta,c,kernel_name):
    # Initialization
    #R = 0.006 # Tensile instability term
    kernel_name = globals()['%s' % kernel_name]
    pressure = array([0.,0.,0.])
    direction = array([1.,1.,1.])
    count = 0

    neighbors = { key: dict_moving[key] for key in dict_moving[i]['Fluid Neighbors'] }
    for j in dict_moving[i]['Fluid Neighbors']:
        r_vector = array(dict_moving[i]['moving_r'][count])
        r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)
        Grad_W = kernel_name(r,h).Gradient()
        direction = sign(around(r_vector, decimals = 10))
        PI = Artificial_Viscosity(dict_moving[i],neighbors[j],alpha,beta,h,c,r_vector)
        pressure = pressure - direction * neighbors[j]['Mass']* \
        ((PI + dict_moving[i]['Pressure']/dict_moving[i]['Density']**2)+(neighbors[j]['Pressure']/neighbors[j]['Density']**2))*Grad_W
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

        direction = sign(around(r_vector, decimals = 10))
        
        surface = surface - direction * delta * Laplacian_c * n/abs(n)
        count += 1

    return surface

def Artificial_Viscosity(current,neighbor,alpha,beta,h,c,r_vector):
    v_ab = array([-neighbor['X Velocity']+current['X Velocity'],-neighbor['Y Velocity']+current['Y Velocity'],-neighbor['Z Velocity']+current['Z Velocity']])
    if dot(v_ab,r_vector) < 0:

        c_ab = (c[current['Type']]+c[neighbor['Type']])/2
        rho_ab = (current['Density']+neighbor['Density'])/2

        mu_ab = h*dot(v_ab,r_vector)/(dot(r_vector,r_vector)+0.01*h**2)

        return -alpha*c_ab*mu_ab + beta*mu_ab**2/rho_ab

    else:

        return 0

# Boundary-Fluid Pressure Force
def Boundary_Fluid_Pressure(dict_moving,dict_boundary,h,count,kernel_name):
    kernel_name = globals()['%s' % kernel_name]
    direction = array([0.,0.,0.])
    r_vector = array(dict_moving['boundary_r'][count])
    r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)
    Grad_W = kernel_name(r,h).Gradient()
    direction = sign(around(r_vector, decimals = 10))

    return direction * Grad_W * -dict_moving['Mass']*dict_boundary['psi']*dict_moving['Pressure']/dict_moving['Density']**2

# Boundary-Fluid Friction Force
def Boundary_Fluid_Friction(dict_moving,dict_boundary,h,delta,c,count,kernel_name):
    kernel_name = globals()['%s' % kernel_name]

    v_ab = array([-dict_boundary['X Velocity']+dict_moving['X Velocity'],-dict_boundary['Y Velocity']+dict_moving['Y Velocity'],-dict_boundary['Z Velocity']+dict_moving['Z Velocity']])
    r_vector = array(dict_moving['boundary_r'][count])
    if dot(v_ab,r_vector) < 0:
        direction = array([0.,0.,0.])
        r = sqrt(r_vector[0]**2+r_vector[1]**2+r_vector[2]**2)
        Grad_W = kernel_name(r,h).Gradient()
        for dir_num in range(0,len(direction)):
            if round(r_vector[dir_num],10) == 0:
                direction[dir_num] = 0
            else: 
                direction[dir_num] = round(r_vector[dir_num],10)/abs(round(r_vector[dir_num],10))

        c_ab = (c[dict_moving['Type']]+c[dict_boundary['Type']])/2

        mu_ab = delta*h*c_ab/(2*dict_moving['Density'])

        PI = -mu_ab*dot(v_ab,r_vector)/(dot(r_vector,r_vector)+0.01*h**2)

        return -dict_moving['Mass']*dict_boundary['psi']*PI*Grad_W*direction
    else:

        return 0

