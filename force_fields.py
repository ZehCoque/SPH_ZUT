import kernels

# Smoothed density
def Density(mass_j,velocity_i,velocity_j,Grad_W):
    density = 0
    for i in range(0,len(mass_j)):
        density = density + mass_j[i]*(velocity_i-velocity_j[i])*Grad_W
    
    return density

# Smoothed pressure
# def Pressure():
#     pressure = 0.
#     for i in range(0,mass_j.__len__):
#         p_tmp = mass_j[i]*(pressure_i[i]+pressure_j[i])/(2*density_j[i])
#         W = kernels.Spiky(posi,posj[i],h)
#         pressure = pressure + p_tmp * W
    
#     return pressure