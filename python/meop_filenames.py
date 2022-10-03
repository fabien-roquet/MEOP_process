from pathlib import Path
import os
import shutil



version = 'MEOP-CTD_2021-11-26'


processdir = Path.home() / 'MEOP_process'
datadir = Path.home() / 'MEOP_dropbox' / 'RAW_MEOP_DATA/'
processdir = Path('/media/disk2/roquet') / 'MEOP_process'
matlabdir = processdir / 'matlab'
inputdir = Path('/home/smru/datadir/all/')
refdir = Path.home() / 'MEOP_dropbox' / 'REF_DATASETS'
publicdir = Path.home() / 'MEOP_process' / 'public'
publicdir_CTD = publicdir / version

# 1. utils to reconstruct name of ncfiles
def deployment_from_smru_name(smru_name):
    return smru_name.split("-")[0]

def smru_name_from_fname_prof(fname_prof):
    return fname_prof.name.split('_')[0]

# return ncARGO filename corresponding to a smru_name
def fname_prof(smru_name,deployment='',qf='lr0'):
    if not deployment:
        deployment = deployment_from_smru_name(smru_name)
    return Path(processdir,'final_dataset_prof',deployment,smru_name+'_'+qf+'_prof.nc')

# return ncARGO filename corresponding to a smru_name
def fname_traj(smru_name,deployment=''):
    if not deployment:
        deployment = deployment_from_smru_name(smru_name)
    return Path(processdir,'final_dataset_prof',deployment,smru_name+'_traj.nc')

# return list of ncARGO filename
def list_fname_prof(smru_name='',deployment='',qf='*',folder=processdir):
    if not deployment:
        deployment = deployment_from_smru_name(smru_name)
    dirEXP = Path(folder) / 'final_dataset_prof' / deployment
    if smru_name:
        prefix = smru_name
    else:
        prefix = deployment+'-*'
    list_fname = [ncfile for ncfile in dirEXP.glob(f'{prefix}_{qf}_prof.nc')]
    return list(set(list_fname))

# return smru_name
def list_smru_name(smru_name='',deployment='',qf='*'):
    list_smru_name = [ncfile.name.split('_')[0] for ncfile in list_fname_prof(smru_name,deployment,qf)]
    return list(set(list_smru_name))

# return ncARGO filename
def fname_plots(smru_name,deployment='',qf='lr0',suffix='_plot'):
    if not deployment:
        deployment = deployment_from_smru_name(smru_name)
    return processdir / 'plots' / deployment / (smru_name+'_'+qf+'_'+suffix+'.png')

# return ncARGO filename
def list_fname_plots(smru_name='',deployment='',qf='*',suffix='_plot'):
    if not deployment:
        deployment = deployment_from_smru_name(smru_name)
    dirEXP = processdir / 'plots' / deployment
    if smru_name:
        prefix = smru_name
    else:
        prefix = deployment+'-*'
    list_fname = [ncfile for ncfile in dirEXP.glob(f'{prefix}_{qf}_{suffix}.png')]
    return list_fname

# copy a file
def copy_file(file_name,src_dir,dst_dir):
    shutil.copyfile(Path(src_dir)/file_name,Path(dst_dir)/file_name)




