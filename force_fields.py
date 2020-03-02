from kernels import *
from numpy import array, sqrt, dot, sign, around
import hashing

# Smoothed density
def Density(current,neighbor,r,h,r_vector=[],beta=[],correction=False,kernel_name='Poly_6'):
    # Initialization
    kernel_name = globals()['%s' % kernel_name]

    if correction == True:
        W = (beta[0]+beta[1]*r_vector[0]+beta[2]*r_vector[1]+beta[3]*r_vector[2])*kernel_name(r,h).Kernel()
    else:
        W = kernel_name(r,h).Kernel()

        return neighbor['Mass']*W

# Smoothed pressure force
def Pressure(current,neighbor,r,h,t,r_vector,kernel_name,alpha=0,beta=0,c=0):
    # Initialization
    #R = 0.006 # Tensile instability term
    kernel_name = globals()['%s' % kernel_name]

    Grad_W = kernel_name(r,h).Gradient(r_vector)
    # PI = Artificial_Viscosity(current,neighbor,alpha,beta,h,c,r_vector)
    # return -neighbor['Mass']*current['Density']* \
    # ((current['Pressure']/current['Density']**2)+(neighbor['Pressure']/neighbor['Density']**2))*Grad_W
    # return array(-current['Mass']/current['Density']*(current['Pressure']+neighbor['Pressure'])/2*Grad_W)

    return array(-current['Mass']/current['Density']*neighbor['Mass']/neighbor['Density']*(current['Pressure']+neighbor['Pressure'])/2*Grad_W) # Force -> must be divided by mass
    # return array(-current['Mass']*neighbor['Mass']*(current['Pressure']/current['Density']**2+neighbor['Pressure']/neighbor['Density']**2)*Grad_W)

# Smoothed Viscosity

def Viscosity(current,neighbor,r,h,t,kernel_name,mu):
    #Initialization
    kernel_name = globals()['%s' % kernel_name]

    Laplacian_W = kernel_name(r,h).Laplacian()
    dv = array([neighbor['X Velocity']-current['X Velocity'],neighbor['Y Velocity']-current['Y Velocity'],neighbor['Z Velocity']-current['Z Velocity']])

    # return array(mu*dv*neighbor['Mass']/neighbor['Density']*Laplacian_W)

    return array(current['Mass']/current['Density']*mu*neighbor['Mass']/neighbor['Density']*dv*Laplacian_W) # Force -> must be divided by mass
        

# Smoothed Surface Tension Force
def Surface_Tension(current,neighbor,r,h,r_vector,kernel_name,delta,rho_0):
    # Initialization
    kernel_name = globals()['%s' % kernel_name]
    W = kernel_name(r,h).Kernel()
    
    # Cohesion force
    cohesion = -delta*current['Mass']*neighbor['Mass']*W*r_vector/r

    # Surface area minimization
    curvature = -delta*current['Mass']*(current['Normal']-neighbor['Normal'])

    # correction factor
    K = 2*rho_0/(current['Density']+neighbor['Density'])

    # Surface tension force
    return K*(cohesion+curvature)
    # if sqrt(dot(current['Normal'],current['Normal'])) > threshold:
    #     kappa = current['Color Map Laplacian']/sqrt(dot(current['Normal'],current['Normal']))
    #     return array(-delta*kappa*neighbor['Mass']*r_vector*W) # Force -> must be divided by mass
    # else:
    #     return array([0.,0.,0.])

    # return array(-delta * n/sqrt(dot(n,n))*neighbor['Mass']/neighbor['Density'] * Laplacian_W)
    
# def Artificial_Viscosity(current,neighbor,alpha,beta,h,c,r_vector):
#     v_ab = array([-neighbor['X Velocity']+current['X Velocity'],-neighbor['Y Velocity']+current['Y Velocity'],-neighbor['Z Velocity']+current['Z Velocity']])
#     if dot(v_ab,r_vector) < 0:

#         c_ab = (c[current['Type']]+c[neighbor['Type']])/2
#         rho_ab = (current['Density']+neighbor['Density'])/2

#         mu_ab = h*dot(v_ab,r_vector)/(dot(r_vector,r_vector)+0.01*h**2)

#         return -alpha*c_ab*mu_ab + beta*mu_ab**2/rho_ab

#     else:

#         return 0

# Boundary-Fluid Pressure Force
def Boundary_Fluid_Pressure(current,neighbor,r,h,r_vector,kernel_name):
    kernel_name = globals()['%s' % kernel_name]
    Grad_W = kernel_name(r,h).Gradient(r_vector)
    #array(Grad_W * -neighbor['psi']*current['Pressure']/current['Density']**2)
    return array(-current['Density']*(current['Pressure']/(current['Density'])**2)*neighbor['Mass']*Grad_W)
    # return array(-Grad_W * neighbor['psi']*current['Pressure']/current['Density']**2)

# Boundary-Fluid Friction Force
def Boundary_Fluid_Friction(current,neighbor,r,h,r_vector,kernel_name,delta,c):
    kernel_name = globals()['%s' % kernel_name]

    v_ab = array([current['X Velocity']-neighbor['X Velocity'],current['Y Velocity']-neighbor['Y Velocity'],current['Z Velocity']-neighbor['Z Velocity']])

    Grad_W = kernel_name(r,h).Gradient(r_vector)

    c_ab = (c[current['Type']]+c[neighbor['Type']])/2

    ni_ab = delta*h*c_ab/(2*current['Density'])

    PI = -ni_ab*max(0,dot(v_ab,r_vector))/(r**2+0.01*h**2)

    return array(-current['Density']*current['Mass']*neighbor['Mass']*PI*Grad_W)

