from pyevtk.hl import pointsToVTK, VtkGroup, VtkData
from pathlib import Path
from os import listdir
from numpy import sqrt, asarray,array, zeros
import csv
import time
import pandas

def get_paths(main_folder):
    '''Make directories and return the paths for the vtk, group and csv folders'''

    Path(main_folder).mkdir(parents=True, exist_ok=True)
    path = main_folder

    folders = listdir(path)

    num = 1

    if not folders:
        foldername = 'simulation 1/'
    else:
        for string in folders:
            if [int(s) for s in string.split() if s.isdigit()][0] > num:
                num = [int(s) for s in string.split() if s.isdigit()][0]
            
        foldername = 'simulation ' + str(num+1) + '/'

    Path(path + foldername + '/csv').mkdir(parents=True, exist_ok=True)
    Path(path + foldername + '/vtk').mkdir(parents=True, exist_ok=True)
    group = VtkGroup(path + foldername + 'moving_particles')

    return [path + foldername + '/vtk',group,path + foldername + '/csv']

def save_csv(path,iteration,dictionary):
    filename = path + '/iter_' + str(iteration) + '.csv' #csv filename
    try:
        headers = dictionary[next(iter(dictionary.keys()))].keys()
    except:
        headers = dictionary.keys()
    try:
        with open(filename, 'w') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=headers,lineterminator = '\n')
            writer.writeheader()
            for i in dictionary:
                writer.writerow(dictionary[i])

    except IOError:
        print()
        print('-'*40)
        print("Error on function save_csv")
        print()

def save_moving_vtk(path,iteration,dictionary):
    filename = path + '/iter_' + str(iteration)  #vtk filename

    X = asarray([dictionary[d].get('X') for d in dictionary])
    Y = asarray([dictionary[d].get('Y') for d in dictionary])
    Z = asarray([dictionary[d].get('Z') for d in dictionary])
    Vx = asarray([dictionary[d].get('X Velocity') for d in dictionary])
    Vy = asarray([dictionary[d].get('Y Velocity') for d in dictionary])
    Vz = asarray([dictionary[d].get('Z Velocity') for d in dictionary])
    V_Mag = array(sqrt(Vx**2+Vy**2+Vz**2))
    # V_vector = asarray([[dictionary[d].get('X Velocity'), dictionary[d].get('Y Velocity'),dictionary[d].get('Z Velocity')] for d in dictionary])
    rho = asarray([dictionary[d].get('Density') for d in dictionary])
    p = asarray([dictionary[d].get('Pressure') for d in dictionary])
    try:
        # Total forces
        total_force_x = asarray([dictionary[d].get('Total Force')[0] for d in dictionary])
        total_force_y = asarray([dictionary[d].get('Total Force')[1] for d in dictionary])
        total_force_z = asarray([dictionary[d].get('Total Force')[2] for d in dictionary])
        # Pressure forces
        p_force_x = asarray([dictionary[d].get('Pressure Force')[0] for d in dictionary])
        p_force_y = asarray([dictionary[d].get('Pressure Force')[1] for d in dictionary])
        p_force_z = asarray([dictionary[d].get('Pressure Force')[2] for d in dictionary])
        # # Boundary-Fluid Pressure forces
        # fb_p_force_x = asarray([dictionary[d].get('Boundary-Fluid Pressure')[0] for d in dictionary])
        # fb_p_force_y = asarray([dictionary[d].get('Boundary-Fluid Pressure')[1] for d in dictionary])
        # fb_p_force_z = asarray([dictionary[d].get('Boundary-Fluid Pressure')[2] for d in dictionary])
        # # Boundary-Fluid Friction forces
        # fb_f_force_x = asarray([dictionary[d].get('Boundary-Fluid Friction')[0] for d in dictionary])
        # fb_f_force_y = asarray([dictionary[d].get('Boundary-Fluid Friction')[1] for d in dictionary])
        # fb_f_force_z = asarray([dictionary[d].get('Boundary-Fluid Friction')[2] for d in dictionary])
        # Viscosity forces
        v_force_x = asarray([dictionary[d].get('Viscosity Force')[0] for d in dictionary])
        v_force_y = asarray([dictionary[d].get('Viscosity Force')[1] for d in dictionary])
        v_force_z = asarray([dictionary[d].get('Viscosity Force')[2] for d in dictionary])
        # Surface tension forces
        st_force_x = asarray([dictionary[d].get('Surface Tension Force')[0] for d in dictionary])
        st_force_y = asarray([dictionary[d].get('Surface Tension Force')[1] for d in dictionary])
        st_force_z = asarray([dictionary[d].get('Surface Tension Force')[2] for d in dictionary])

    except:
        # Total forces
        total_force_x = zeros(len(X))
        total_force_y = zeros(len(X))
        total_force_z = zeros(len(X))
        # Pressure forces
        p_force_x = zeros(len(X))
        p_force_y = zeros(len(X))
        p_force_z = zeros(len(X))
        # # Boundary-Fluid Pressure forces
        # fb_p_force_x = zeros(len(X))
        # fb_p_force_y = zeros(len(X))
        # fb_p_force_z = zeros(len(X))
        # # Boundary-Fluid Friction forces
        # fb_f_force_x = zeros(len(X))
        # fb_f_force_y = zeros(len(X))
        # fb_f_force_z = zeros(len(X))
        # Viscosity forces
        v_force_x = zeros(len(X))
        v_force_y = zeros(len(X))
        v_force_z = zeros(len(X))
        # Surface tension forces
        st_force_x = zeros(len(X))
        st_force_y = zeros(len(X))
        st_force_z = zeros(len(X))
        

    pointsToVTK(filename,X,Y,Z, 
    data = {"Vx" : Vx, "Vy" : Vy, "Vz" : Vz,"V_Mag":V_Mag, "rho" : rho,"p" : p,"Total Force X" : total_force_x, "Total Force Y" : total_force_y, "Total Force Z" : total_force_z,
    "Pressure X" : p_force_x,"Pressure Y" : p_force_y,"Pressure Z" : p_force_z,"Viscosity X" : v_force_x,"Viscosity Y" : v_force_y,"Viscosity Z" : v_force_z,
    "ST X" : st_force_x,"ST Y" : st_force_y,"ST Z" : st_force_z})
    
def save_boundary_vtk(path,dictionary):
    filename = path.replace('/vtk','')  + '/boundary' #vtk filename
    pointsToVTK(filename, asarray([dictionary[d].get('X') for d in dictionary]), asarray([dictionary[d].get('Y') for d in dictionary]), asarray([dictionary[d].get('Z') for d in dictionary]),
    data = {"mass" : asarray([dictionary[d].get('Mass') for d in dictionary]) })

def add_to_group(path,iteration,time,group):
    filename = path + '/iter_' + str(iteration)  + '.vtu' #vtk filename 
    group.addFile(filepath = filename, sim_time = time)

def save_group(group):
    group.save()

def info(current_time, final_time, start_time, now,delta_t,iteration,max_rho_err, message='Simulating...'):
        print(
            "\r{} |{}{}| {}% -/-/- Current simulation time: {}s -/- Time step: {}s -/- Elapsed time: {}s -/- Iteration: {} -/- Max Density Error: {}"
            .format(
                message,
                "#" * int(20 * current_time / final_time),
                " " * (20 - int(20 * current_time / final_time)),
                int(100 * current_time / final_time),
                round(current_time,6),
                round(delta_t,6),
                round(abs(start_time - now)),
                iteration,
                round(max_rho_err,3)
            ),
            flush=True,
            end="  "
        )

def csv_to_dict(path):
    df = pandas.read_csv(path)
    mydict = { } 
    for i in range(0,df.shape[0]):
        mydict[i] = { }
    for key in df:
        array = df[key].values
        for j in range(0,len(array)):
            mydict[j].update({key: array[j]})

    return mydict

def continue_last_sim():
    folders = listdir("./results")
    max_num1 = 0
    for i in folders:
        num = [int(s) for s in i.split() if s.isdigit()][0]
        if num > max_num1:
            max_num1 = num
    
    sim_folder = "./results/simulation " + str(max_num1) + "/csv/"
    files = listdir(sim_folder)
    max_num2 = 0
    for i in files:
        num = int(i.split("_")[1].replace(".csv",""))
        if num > max_num2:
            max_num2 = num

    return [csv_to_dict(sim_folder + "iter_" + str(max_num2 - 1) + ".csv"), max_num2-1, "./results/simulation " + str(max_num1), VtkGroup("./results/simulation " + str(max_num1) + "/moving_particles")]