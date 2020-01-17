from pyevtk.hl import pointsToVTK, VtkGroup
from pathlib import Path
from os import listdir

def get_paths(main_folder):
    '''Make directories and return the paths for the vtk, group and csv folders'''

    Path("./results/").mkdir(parents=True, exist_ok=True)
    path = "./results/"

    folders = listdir(path)

    if not folders:
        foldername = 'simulation 1/'
    else:
        for string in folders:
            num = int(string[-1])
            
        foldername = 'simulation ' + str(num+1) + '/'

    Path(path + foldername + '/csv').mkdir(parents=True, exist_ok=True)
    Path(path + foldername + '/vtk').mkdir(parents=True, exist_ok=True)
    group = VtkGroup(path + foldername + 'vtk')

    return [path + foldername + '/vtk',group,path + foldername + '/csv']

def save_csv(path,iter_time,df):
    filename = path + '/sph_' + str(round(iter_time,3)) + '.csv' #csv filename
    df.to_csv(filename,index=False) #csv save

def save_vtk(path,iter_time,df):
    filename = path + '/sph_' + str(round(iter_time,3)) #vtk filename
    pointsToVTK(filename, df['X'].values, df['Y'].values, df['Z'].values, 
    data = {"Vx" : df['X Velocity'].values, "Vy" : df['Y Velocity'].values, "Vz" : df['Z Velocity'].values,
    "rho" : df['Density'].values,"p" : df['Pressure'].values})

def add_to_group(path,time,group):
    filename = path + '/sph_' + str(round(time,3)) + '.vtu' #vtk filename 
    group.addFile(filepath = filename, sim_time = time)

def save_group(group):
    group.save()

