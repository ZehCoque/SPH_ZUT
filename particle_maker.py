from numpy import linspace, array
import sys

def make_prism(x,y,z,num_x,num_y,num_z,radius,mass,type_,vx=0.,vy=0.,vz=0.,dict_index=0,prism={}):
    ''' returns a prism given the initial coordinates and the number of particles in each axis \n
    x, y and z are the coordinates of the first particle of the prism. \n
    num_x, num_y and num_z are the numbers of particles for each axis (first particle included).
    They can be either positive or negative. Negative numbers plot particles in the oposite direction of the global axis. \n
    radius, mass and type stands for the radius, mass and type for all particles in the prism. \n
    vx, vy and vz are the velocity coordinates for all the particles in the prism. This is optional and the dafult value is 0. \n
    dict_index is the first index of the dictionary if one wants to continue to build from a previous dictionary \n
    prism is a previous dictionary'''

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
                prism[dict_index] = {'X': i,
                                     'Y': j,
                                     'Z': k,
                                     'X Velocity':vx,
                                     'Y Velocity':vy,
                                     'Z Velocity':vz,
                                     'Pressure': 0.,
                                     'Density': 0.,
                                     'Mass':mass,
                                     'Type':type_}
                dict_index = dict_index + 1

    return prism

def make_prism2(coord_init,coord_final,radius,mass,type_,vx=0.,vy=0.,vz=0.,dict_index=0,prism={}):
    ''' returns a prism given the initial and final coordinates \n
    coord_init -> array(x,y,z) of the initial coordinates \n
    coord_final -> array(x,y,z) of the final coordinates \n
    See make_prism function for the other variables'''
    
    
    distance = array(coord_final) - array(coord_init)

    if distance[0]/(2*radius)%1 != 0 and distance[1]/(2*radius)%1 != 0 and distance[2]/(2*radius)%1 != 0:
        print()
        print("-" * 40)
        print('Division by radius did not return a whole number while trying to make a prism')    
        sys.exit(1)
    else:
        
        num_x = int((coord_final[0] - coord_init[0])/(2*radius))
        num_y = int((coord_final[1] - coord_init[1])/(2*radius))
        num_z = int((coord_final[2] - coord_init[2])/(2*radius))

        return make_prism(coord_init[0],coord_init[1],coord_init[2],num_x,num_y,num_z,radius,mass,type_,vx,vy,vz,dict_index,prism)

def make_box(coord_init,coord_final,radius,mass,type_,vx=0.,vy=0.,vz=0.,dict_index=0,box={}):
    ''' returns an empty box given the initial and final coordinates \n
    coord_init -> array(x,y,z) of the initial coordinates \n
    coord_final -> array(x,y,z) of the final coordinates across main diagonal \n
    See make_prism function for the other variables'''

    distance = array(coord_final) - array(coord_init)

    if distance[0]/(2*radius)%1 != 0 and distance[1]/(2*radius)%1 != 0 and distance[2]/(2*radius)%1 != 0:
        print()
        print("-" * 40)
        print('Division by radius did not return a whole number while trying to make a prism')    
        sys.exit(1)
    else:
        
        num_x = int((coord_final[0] - coord_init[0])/(2*radius))
        num_y = int((coord_final[1] - coord_init[1])/(2*radius))
        num_z = int((coord_final[2] - coord_init[2])/(2*radius))

        vx = float(vx)
        vy = float(vy)
        vz = float(vz)    

        x_array = [float(coord_init[0])]
        for i in range(1,abs(num_x)):
            x_array.append(float(x_array[-1]+2*radius))
        y_array = [float(coord_init[1])]
        for i in range(1,abs(num_y)):
            y_array.append(float(y_array[-1]+2*radius))
        z_array = [float(coord_init[2])]
        for i in range(1,abs(num_z)):
            z_array.append(float(z_array[-1]+2*radius))

        for i in x_array:
            for j in y_array:
                for k in z_array:
                    if round(i,3) in coord_init or round(j,3) in coord_init or round(k,3) in coord_init:
                        box[dict_index] = {'X': i,
                                           'Y': j,
                                           'Z': k,
                                           'X Velocity':vx,
                                           'Y Velocity':vy,
                                           'Z Velocity':vz,
                                           'Pressure': 0.,
                                           'Density': 0.,
                                           'Mass':mass,
                                           'Type':type_}
                        dict_index = dict_index + 1

        for i in x_array:
            for j in y_array:
                for k in z_array:
                    if (round(i+2*radius,3) in coord_final or round(j+2*radius,3) in coord_final or round(k+2*radius,3) in coord_final) and \
                    (round(i,3) not in coord_init and round(j,3) not in coord_init and round(k,3) not in coord_init):
                        box[dict_index] = {'X': i,
                                           'Y': j,
                                           'Z': k,
                                           'X Velocity':vx,
                                           'Y Velocity':vy,
                                           'Z Velocity':vz,
                                           'Pressure': 0.,
                                           'Density': 0.,
                                           'Mass':mass,
                                           'Type':type_}
                        dict_index = dict_index + 1

        return box