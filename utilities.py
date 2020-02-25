from pyevtk.hl import pointsToVTK, VtkGroup
from pathlib import Path
from os import listdir
from numpy import sqrt, asarray
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
                try:
                    writer.writerow(dictionary[i])
                except:
                    writer.writerow(dictionary)
                    break

    except IOError:
        print()
        print('-'*40)
        print("Error on function save_csv")
        print()

def save_moving_vtk(path,iteration,dictionary):
    filename = path + '/iter_' + str(iteration)  #vtk filename
    # try:
    pointsToVTK(filename, asarray([dictionary[d].get('X') for d in dictionary]), asarray([dictionary[d].get('Y') for d in dictionary]), asarray([dictionary[d].get('Z') for d in dictionary]), 
    data = {"Vx" : asarray([dictionary[d].get('X Velocity') for d in dictionary]), "Vy" : asarray([dictionary[d].get('Z Velocity') for d in dictionary]), "Vz" : asarray([dictionary[d].get('Z Velocity') for d in dictionary]),
    "rho" : asarray([dictionary[d].get('Density') for d in dictionary]),"p" : asarray([dictionary[d].get('Pressure') for d in dictionary])})
    
def save_boundary_vtk(path,dictionary):
    filename = path.replace('/vtk','')  + '/boundary' #vtk filename
    pointsToVTK(filename, asarray([dictionary[d].get('X') for d in dictionary]), asarray([dictionary[d].get('Y') for d in dictionary]), asarray([dictionary[d].get('Z') for d in dictionary]),
    data = {"psi" : asarray([dictionary[d].get('psi') for d in dictionary]) })

def add_to_group(path,iteration,time,group):
    filename = path + '/iter_' + str(iteration)  + '.vtu' #vtk filename 
    group.addFile(filepath = filename, sim_time = time)

def save_group(group):
    group.save()

def info(current_time, final_time, start_time, now,delta_t,iteration, message='Simulating...'):
        print(
            "\r{} |{}{}| {}% -/-/- Current simulation time: {}s -/-/- Time step: {}s -/-/- Elapsed time: {}s -/-/- Iteration: {}"
            .format(
                message,
                "#" * int(20 * current_time / final_time),
                " " * (20 - int(20 * current_time / final_time)),
                int(100 * current_time / final_time),
                round(current_time,6),
                round(delta_t,6),
                round(abs(start_time - now)),
                iteration
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