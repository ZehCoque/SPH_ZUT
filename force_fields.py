from kernels import *
from numpy import array, sqrt, dot, sign, around
import hashing

# Smoothed density
def Density(current,neighbor,r,h,t,r_vector=[],beta=[],correction=False,kernel_name='Poly_6'):
    # Initialization
    kernel_name = globals()['%s' % kernel_name]

    if correction == True:
        W = (beta[0]+beta[1]*r_vector[0]+beta[2]*r_vector[1]+beta[3]*r_vector[2])*kernel_name(r,h).Kernel()
    else:
        W = kernel_name(r,h).Kernel()

    if t == 1: # Fluid - Fluid Interactions
        return neighbor['Mass']*W
    elif t == 0: # Fluid - Boundary Interactions
        return neighbor['psi']*W
    else:
        return "error"

# Smoothed pressure force
def Pressure(current,neighbor,r,h,t,r_vector,kernel_name,alpha=0,beta=0,c=0):
    # Initialization
    #R = 0.006 # Tensile instability term
    kernel_name = globals()['%s' % kernel_name]
    pressure = array([0.,0.,0.])
    direction = sign(around(r_vector, decimals = 10))

    Grad_W = kernel_name(r,h).Gradient()
    # PI = Artificial_Viscosity(current,neighbor,alpha,beta,h,c,r_vector)
    pressure = -direction *  neighbor['Mass']* \
    ((current['Pressure']/current['Density']**2)+(neighbor['Pressure']/neighbor['Density']**2))*Grad_W

    return array(pressure)

# Smoothed Viscosity

def Viscosity_Kernel(current,neighbor,r,h,t,kernel_name,mu):
    #Initialization
    kernel_name = globals()['%s' % kernel_name]
    viscosity = array([0.,0.,0.])

    Laplacian_W = kernel_name(r,h).Laplacian()
    dv = array([neighbor['X Velocity']-current['X Velocity'],neighbor['Y Velocity']-current['Y Velocity'],neighbor['Z Velocity']-current['Z Velocity']])

    viscosity = neighbor['Mass']*dv/neighbor['Density']*Laplacian_W
        
    return array(mu*viscosity)

# Smoothed Surface Tension Force
def Surface_Tension(current,neighbor,r,h,t,r_vector,kernel_name,delta):
    # Initialization
    threshold =  7.065
    kernel_name = globals()['%s' % kernel_name]
    direction = sign(around(r_vector, decimals = 10))

    Grad_W = kernel_name(r,h).Gradient()
    n = neighbor['Mass']/neighbor['Density'] * Grad_W * direction

    if sqrt(n[0]**2+n[1]**2+n[2]**2) >= threshold:
    
        Laplacian_W = kernel_name(r,h).Laplacian()
        Laplacian_c = neighbor['Mass']/neighbor['Density'] * Laplacian_W
        
        return array(delta * Laplacian_c * n/(sqrt(n[0]**2+n[1]**2+n[2]**2)))
    
    else:
        return 0

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
def Boundary_Fluid_Pressure(current,neighbor,r,h,r_vector,kernel_name):
    kernel_name = globals()['%s' % kernel_name]
    direction = sign(around(r_vector, decimals = 10))
    Grad_W = kernel_name(r,h).Gradient()
    direction = sign(around(r_vector, decimals = 10))

    return array(direction * Grad_W * -current['Mass']*neighbor['psi']*current['Pressure']/current['Density']**2)

# Boundary-Fluid Friction Force
def Boundary_Fluid_Friction(current,neighbor,r,h,r_vector,kernel_name,delta,c):
    kernel_name = globals()['%s' % kernel_name]

    v_ab = array([-neighbor['X Velocity']+current['X Velocity'],-neighbor['Y Velocity']+current['Y Velocity'],-neighbor['Z Velocity']+current['Z Velocity']])

    if dot(v_ab,r_vector) < 0:

        Grad_W = kernel_name(r,h).Gradient()
        direction = sign(around(r_vector, decimals = 10))

        c_ab = (c[current['Type']]+c[neighbor['Type']])/2

        mu_ab = delta*h*c_ab/(2*current['Density'])

        PI = -mu_ab*dot(v_ab,r_vector)/(abs(r_vector)**2+0.01*h**2)

        return array(-current['Mass']*neighbor['psi']*PI*Grad_W*direction)
    else:

        return array([0.,0.,0.])

