from numpy import linspace

def make_prism(x,y,z,num_x,num_y,num_z,radius,mass,type_,vx=0.,vy=0.,vz=0.,count=0,prism={}):
    ''' filaname is the name of the file this function will save the particles. \n
    x, y and z are the coordinates of the first particle of the prism. \n
    num_x, num_y and num_z are the numbers of particles for each axis (first particle included).
    They can be either positive or negative. Negative numbers plot particles in the oposite direction of the global axis. \n
    radius, mass and type stands for the radius, mass and type for all particles in the prism. \n
    vx, vy and vz are the velocity coordinates for all the particles in the prism. This is optional and the dafult value is 0. \n'''

    vx = float(vx)
    vy = float(vy)
    vz = float(vz)

    x_array = [float(x)]
    for i in range(1,abs(num_x)):
        x_array.append(float(x_array[-1]+2*radius))
    y_array = [float(y)]
    for i in range(1,abs(num_y)):
        y_array.append(float(y_array[-1]+2*radius))
    z_array = [float(z)]
    for i in range(1,abs(num_z)):
        z_array.append(float(z_array[-1]+2*radius))

    for i in x_array:
        for j in y_array:
            for k in z_array:
                prism[count] = {'X': i,
                'Y': j,
                'Z': k,
                'X Velocity':vx,
                'Y Velocity':vy,
                'Z Velocity':vz,
                'Pressure': 0.,
                'Density': 0.,
                'Mass':mass,
                'Type':type_}
                count = count + 1

    return prism
