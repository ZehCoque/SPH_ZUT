import kernels

# Smoothed pressure force field
def Pressure(posi=[0,0,0],posj=[[0,0,0]],h=2,mass_j=[0],pressure_i=0,pressure_j=[0],density_j=[0]):
    pressure = 0.
    for i in range(0,mass_j.__len__):
        p_tmp = mass_j[i]*(pressure_i[i]+pressure_j[i])/(2*density_j[i])
        W = kernels.Spiky(posi,posj[i],h)
        pressure = pressure + p_tmp * W
    
    return pressure    