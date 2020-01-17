from pyevtk.hl import pointsToVTK, VtkGroup
from pathlib import Path
from os import listdir
from numpy import sqrt

def get_paths(main_folder):
    '''Make directories and return the paths for the vtk, group and csv folders'''

    Path(main_folder).mkdir(parents=True, exist_ok=True)
    path = main_folder

    folders = listdir(path)

    if not folders:
        foldername = 'simulation 1/'
    else:
        for string in folders:
            try:
                num = int(string[-1])
            except:
                num = num
            
        foldername = 'simulation ' + str(num+1) + '/'

    Path(path + foldername + '/csv2').mkdir(parents=True, exist_ok=True) # temporary
    Path(path + foldername + '/csv').mkdir(parents=True, exist_ok=True)
    Path(path + foldername + '/vtk').mkdir(parents=True, exist_ok=True)
    group = VtkGroup(path + foldername + 'vtk')

    return [path + foldername + '/vtk',group,path + foldername + '/csv']

def save_csv(path,iteration,df):
    filename = path + '/iter_' + str(iteration) + '.csv' #csv filename
    df.to_csv(filename,index=False) #csv save

def save_vtk(path,iteration,df):
    filename = path + '/iter_' + str(iteration)  #vtk filename
    pointsToVTK(filename, df['X'].values, df['Y'].values, df['Z'].values, 
    data = {"Vx" : df['X Velocity'].values, "Vy" : df['Y Velocity'].values, "Vz" : df['Z Velocity'].values,
    "V" : sqrt( df['X Velocity'].values**2+ df['Y Velocity'].values**2+ df['Z Velocity'].values**2),
    "rho" : df['Density'].values,"p" : df['Pressure'].values})

def add_to_group(path,iteration,time,group):
    filename = path + '/iter_' + str(iteration)  + '.vtu' #vtk filename 
    group.addFile(filepath = filename, sim_time = time)

def save_group(group):
    group.save()

