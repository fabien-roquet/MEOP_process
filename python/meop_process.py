# matlab.engine must be imported first
import matlab.engine

from pathlib import Path
import os
import shutil
import xarray as xr
import pandas as pd
import numpy as np
import gsw
import netCDF4
from importlib import reload
import io
import meop


processdir = meop.processdir

#--------------------  MATLAB  ----------------------------#

# pointer to matlab engine
eng = []

def start_matlab():
    # Start matlab iif not already started
    global eng
    try:
        eng.eval("disp('matlab already started')",nargout=0)
        print('matlab already started')
        print('PWD:',eng.pwd())
    except:
        eng = matlab.engine.start_matlab()
        print('matlab started')
        print('PWD:',eng.pwd())
    return

def stop_matlab():
    global eng
    eng.quit()
    eng = []
    return
    

def print_matlab(namevar):
    with io.StringIO() as out:
        eng.eval(namevar,nargout=0,stdout=out,stderr=out)
        print(out.getvalue())
    
def run_command(cmd,verbose=True):
    # execute a matlab command
    with io.StringIO() as out, io.StringIO() as err:
        try:
            print(cmd)
            eng.eval(cmd,nargout=0,stdout=out,stderr=err)
            if verbose:
                print(out.getvalue())
            return True
        except:
            if verbose:
                print(err.getvalue())
            return False
        
# init mirounga and load conf
def init_mirounga():
    eng.addpath(str(processdir))
    conf = eng.eval("init_config();",nargout=1)
    run_command("conf = init_mirounga;")
    return conf


def load_info_deployment(deployment='',smru_name=''):
    init_mirounga()
    eng.workspace['EXP'] = deployment
    eng.workspace['one_smru_name'] = smru_name
    eng.eval("info_deployment=load_info_deployment(conf,EXP,one_smru_name);",nargout=0)

    
def process_tags(deployment='',smru_name=''):
    load_info_deployment(deployment=deployment,smru_name=smru_name)
    if eng.eval("isfield(info_deployment,'invalid_code')") and eng.eval("info_deployment.invalid_code"):
        return False
    if not run_command("remove_deployment(conf,EXP,one_smru_name);"):
        return False
    if not run_command("create_ncargo(conf,EXP,one_smru_name);"):
        return False
    if not run_command("create_fr0(conf,EXP,one_smru_name);"):
        return False
    if not run_command("update_metadata(conf,EXP,one_smru_name);"):
        return False
    if not run_command("apply_adjustments(conf,EXP,one_smru_name);"):
        return False
    if not run_command("apply_tlc(conf,EXP,one_smru_name);"):
        return False
    if not run_command("apply_tlc_fr(conf,EXP,one_smru_name);"):
        return False
    return True


def generate_calibration_plots(deployment='',smru_name=''):
    load_info_deployment(deployment=deployment,smru_name=smru_name)
    if eng.eval("isfield(info_deployment,'invalid_code')") and eng.eval("info_deployment.invalid_code"):
        return False
    if not run_command("generate_plot1(conf,EXP,one_smru_name);"):
        return False
    return True



# update metadata
def update_metadata(deployment='',smru_name=''):

    path_meta = processdir / 'table_meta.csv'
    if not path_meta.exists():
        print(f'Warning: File {path_meta} not found. Metadata not updated.')
        return
    
    df_meta = pd.read_csv(path_meta).set_index('smru_platform_code')

    if smru_name in df_meta.index:
        df_meta = df_meta.loc[[smru_name],:]
    elif deployment:
        df_meta['deployment'] = df_meta.index
        df_meta['deployment'] = df_meta.deployment.str.split('-').apply(lambda x: x[0])
        df_meta = df_meta.loc[df_meta.deployment == deployment,:]

    modes = ['lr0','lr1','hr0','hr1','fr0','fr1']
    for smru_name in df_meta.index:

        meta_row = df_meta.loc[smru_name,:].dropna()

        for qf in modes:

            namefile = meop.fname_prof(smru_name,qf=qf)        
            if Path(namefile).exists():

                with netCDF4.Dataset(namefile,'a') as f:
                    for col in meta_row.keys():
                        if f.location != meta_row[col]:
                            f.location = meta_row[col]

    return


def generate_descriptive_plots(smru_name='',deployment=''):
    import meop
    
    list_qf = ['lr0','hr1','fr1']
    
    for qf in list_qf:

        for smru_name in meop.list_smru_name(smru_name,deployment,qf=qf):
            namefile = meop.fname_prof(smru_name,qf=qf)
            ds = meop.read_ncfile(namefile)
            ds.plot_data_tags('_ADJUSTED',namefig=meop.fname_plots(smru_name,qf=qf,suffix='profiles'))
            ds.plot_TSsections('_ADJUSTED',namefig=meop.fname_plots(smru_name,qf=qf,suffix='sections'))
            ds.close()
    

    
# create a netcdf file combining data on original levels and interp data
def create_final_ncfile(smru_name='',deployment=''):
    
    # copy lr1 in all
    ncfile_lr1 = meop.fname_prof(smru_name,qf='lr1')
    ncfile_hr1 = meop.fname_prof(smru_name,qf='hr1')
    ncfile_out = meop.fname_prof(smru_name,qf='all')

    if ncfile_lr1.is_file() and ncfile_hr1.is_file():
        print('copy: '+ncfile_lr1.name)
        shutil.copyfile(ncfile_hr1,ncfile_out)

    return ncfile_out


# Execute in terminal command line
if __name__ == "__main__":

    import argparse
    
    # create parser
    parser = argparse.ArgumentParser()

    # add arguments to the parser
    parser.add_argument("--smru_name", default ='', help = "Process only SMRU PLATFORM CODE. Value of DEPLOYMENT_CODE not considered.")
    parser.add_argument("--deployment", default ='', help = "Process all tags in DEPLOYMENT_CODE")
    parser.add_argument("--do_all", help = "Process data and produce plots", action='store_true')
    parser.add_argument("--process_data", help = "Process data", action='store_true')
    parser.add_argument("--metadata", help = "Update metadata", action='store_true')
    parser.add_argument("--calibration_plots", help = "Produce calibration plots", action='store_true')
    parser.add_argument("--descriptive_plots", help = "Produce descriptive plots", action='store_true')
    
    # parse the arguments
    args = parser.parse_args()

    smru_name = args.smru_name
    deployment = args.deployment
    
    
    if (smru_name or deployment) and (args.process_data or args.metadata or args.calibration_plots):
        start_matlab()
        conf = init_mirounga()
        if args.process_data or args.do_all:
            process_tags(deployment=deployment,smru_name=smru_name)
        if args.metadata or args.do_all:
            update_metadata(table_meta = table_meta, deployment=deployment,smru_name=smru_name)
        if args.calibration_plots or args.do_all:
            generate_calibration_plots(deployment=deployment,smru_name=smru_name)
        stop_matlab()
    
    if args.descriptive_plots or args.do_all:
        generate_descriptive_plots(smru_name=smru_name,deployment=deployment)
        

        
