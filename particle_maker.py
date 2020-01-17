from numpy import linspace

# def make_prism(filename,x,y,z,num_x,num_y,num_z,radius,mass,type_,vx=0,vy=0,vz=0):
#     ''' filaname is the name of the file this function will save the particles. \n
#     x, y and z are the coordinates of the first particle of the prism. \n
#     num_x, num_y and num_z are the numbers of particles for each axis (first particle included).
#     They can be either positive or negative. Negative numbers plot particles in the oposite direction of the global axis. \n
#     radius, mass and type stands for the radius, mass and type for all particles in the prism. \n
#     vx, vy and vz are the velocity coordinates for all the particles in the prism. This is optional and the dafult value is 0. \n'''

x = 0
y = -2
z = 4
vx = 0
vy = 0
vz = 0
mass = 1
type_ = 'Water'
radius = 0.5
num_x = 2
num_y = 2
num_z = 2

prism = {'X': [],
        'Y': [],
        'Z': [],
            'X Velocity':[],
            'Y Velocity':[],
            'Z Velocity':[],
            'Pressure': [],
            'Density': [],
            'Mass':[],
            'Type':[]}

x_array = [x]
for i in range(1,abs(num_x)):
    x_array.append(x_array[-1]+2*radius)
y_array = [y]
for i in range(1,abs(num_y)):
    y_array.append(y_array[-1]+2*radius)
z_array = [z]
for i in range(1,abs(num_z)):
    z_array.append(z_array[-1]+2*radius)

for i in x_array:
    for j in y_array:
        for k in z_array:
            array = [i,j,k,vx,vy,vz,0,0,mass,type_]
            count = 0
            for var in prism:
                prism[var].append(array[count])
                count = count + 1
                
