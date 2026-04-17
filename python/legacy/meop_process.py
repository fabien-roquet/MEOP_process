# matlab.engine must be imported first
import matlab.engine

from pathlib import Path
import os
import shutil
from shutil import copyfile
import zipfile
import xarray as xr
import pandas as pd
import numpy as np
import gsw
import netCDF4
from importlib import reload
import io
import meop
import meop_filenames
import meop_metadata
_ = os.umask(int('002', 8))

processdir = meop_filenames.processdir


#--------------------  MATLAB  ----------------------------#

# pointer to matlab engine
eng = []

def start_matlab():
    # Start matlab iif not already started
    global eng
    try:
        eng.eval("disp('matlab already started')",nargout=0)
        print('matlab already started')
        run_command(f"cd {processdir};")
        print('PWD:',eng.pwd())
    except:
        eng = matlab.engine.start_matlab()
        print('matlab started')
        run_command(f"cd {processdir};")
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
    eng.addpath(str(processdir / 'matlab'))
    conf = eng.eval("init_config();",nargout=1)
    run_command("conf = init_mirounga;",verbose=False)
    return conf


def load_info_deployment(deployment='',smru_name=''):
    init_mirounga()
    eng.workspace['EXP'] = deployment
    eng.workspace['one_smru_name'] = smru_name
    try :
        eng.eval("info_deployment=load_info_deployment(conf,EXP,one_smru_name);",nargout=0)
    except :
        print("Check file list_deployment.csv for any UNKNOWN country")


def update_config_files(conf):
    
    import filecmp
    from datetime import datetime
    datenow = datetime.now().strftime('%Y%m%d')

    list_files = ['deployment2','deployment3','platform2','platform3']
    
    for file in list_files:
        file_smru  = meop_filenames.inputdir / f'../{file}.json'
        file_local = meop_filenames.datadir / 'config_files' / f'{file}.json'
        file_local_save = meop_filenames.datadir / 'config_files' / f'{file}_{datenow}.json'
        if file_smru.exists():
            if file_local.exists():
                if (not filecmp.cmp(file_smru,file_local)):
                    print(f'{file}.json has been updated')
                    copyfile(file_local,file_local_save)
                    copyfile(file_smru,file_local)
            else:
                print(f'{file}.json has been imported')
                copyfile(file_smru,file_local)
        
    
    return init_mirounga()

    
def import_raw_data(deployment=''):
    
    import filecmp
    
    if not deployment: return
    
    if (meop_filenames.inputdir / deployment).is_dir():
        
        zipfile_orig = meop_filenames.inputdir / deployment / (deployment+'_ODV.zip')
        if not zipfile_orig.exists(): return
        
        zipfile_copy = meop_filenames.datadir / 'raw_smru_data_odv' / (deployment+'_ODV.zip')
        
        if (not zipfile_copy.exists()) or (not filecmp.cmp(zipfile_orig,zipfile_copy)):
            copyfile(zipfile_orig,zipfile_copy)
            with zipfile.ZipFile(zipfile_copy) as z:
                for file in z.namelist():
                    print(file)
                z.extractall(path = meop_filenames.datadir / 'raw_smru_data_odv')

            output = run_command(f'deployment_code = \'{deployment}\';')
            output = run_command('fusion_profilTS_profilFL(deployment_code,conf.rawdir);')
        
    else:
        
        print(f"Folder with raw data is not found: {meop_filenames.inputdir / deployment}")
        
    return


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

    modes = ['lr0','hr0','fr0']
    for smru_name in df_meta.index:
        meta_row = df_meta.loc[smru_name,:].dropna()
        for qf in modes:
            namefile = meop_filenames.fname_prof(smru_name,qf=qf)        
            if Path(namefile).exists():
                with netCDF4.Dataset(namefile,'a') as f:
                    for col in meta_row.keys():
                        if f.location != meta_row[col]:
                            f.location = meta_row[col]

    return


def process_tags(deployment='',smru_name='',notlc=False):
    
    if smru_name:
        print('Process tag :'+smru_name)
        deployment = meop_filenames.deployment_from_smru_name(smru_name)
        print('')
    elif deployment:
        print('Process deployment :'+deployment)
        print('')
    
    load_info_deployment(deployment=deployment,smru_name=smru_name)    
    if eng.eval( "isfield(info_deployment,'invalid_code')") and eng.eval("info_deployment.invalid_code"): return False
    conf = init_mirounga()
    conf = update_config_files(conf)
    import_raw_data(deployment=deployment)
    
    if not run_command("remove_deployment(conf,EXP,one_smru_name);"): return False    
    if not run_command("create_ncargo(conf,EXP,one_smru_name);"): return False    
    if not run_command("create_fr0(conf,EXP,one_smru_name);"): return False    
    if not run_command("create_fr0_without_lr0(conf,EXP,one_smru_name);"): return False    
    if not run_command("update_metadata(conf,EXP,one_smru_name);"): return False    
    update_metadata(deployment=deployment,smru_name=smru_name)    
    if not run_command("apply_adjustments(conf,EXP,one_smru_name);"): return False
    
    if notlc:        
        if not run_command("apply_notlc(conf,EXP,one_smru_name);"): return False
        if not run_command("apply_notlc_fr(conf,EXP,one_smru_name);"): return False
    else:        
        if not run_command("apply_tlc(conf,EXP,one_smru_name);"): return False
        if not run_command("apply_tlc_fr(conf,EXP,one_smru_name);"): return False
    
    return True


def create_hr2(deployment='',smru_name=''):
    
    if smru_name:
        print('Process tag :'+smru_name)
        print('')
    elif deployment:
        print('Process deployment :'+deployment)
        print('')
    
    load_info_deployment(deployment=deployment,smru_name=smru_name)
    if eng.eval("isfield(info_deployment,'invalid_code')") and eng.eval("info_deployment.invalid_code"): return False
    if not run_command("create_hr2(conf,EXP,one_smru_name);"): return False
    return True


def generate_calibration_plots(deployment='',smru_name=''):
    load_info_deployment(deployment=deployment,smru_name=smru_name)
    if eng.eval("isfield(info_deployment,'invalid_code')") and eng.eval("info_deployment.invalid_code"): return False

    if not run_command("generate_plot1(conf,EXP,one_smru_name);"): return False
    return True


def generate_doc_latex(deployment='',smru_name=''):
    load_info_deployment(deployment=deployment,smru_name=smru_name)
    if eng.eval("isfield(info_deployment,'invalid_code')") and eng.eval("info_deployment.invalid_code"): return False

    if not run_command("generate_plot2(conf,EXP,one_smru_name);"): return False
    return True

def export_odv4(deployment='',smru_name=''):
    load_info_deployment(deployment=deployment,smru_name=smru_name)
    if eng.eval("isfield(info_deployment,'invalid_code')") and eng.eval("info_deployment.invalid_code"): return False

    if not run_command("generate_odv4(conf,EXP,one_smru_name);"): return False
    return True


# Execute in terminal command line
if __name__ == "__main__":

    import argparse
    
    # create parser
    parser = argparse.ArgumentParser()

    # add arguments to the parser
    parser.add_argument("--smru_name", default ='', help = "Process only SMRU PLATFORM CODE.")
    parser.add_argument("--deployment", default ='', help = "Process all tags in DEPLOYMENT_CODE")
    parser.add_argument("--do_all", help = "Process data and produce plots", action='store_true')
    parser.add_argument("--update_config_files", help = "Update config json files", action='store_true')
    parser.add_argument("--import_data", help = "Import raw data", action='store_true')
    parser.add_argument("--process_data", help = "Process data", action='store_true')
    parser.add_argument("--create_hr2", help = "Create a netcdf combining hr1 and fr1", action='store_true')
    parser.add_argument("--calibration_plots", help = "Produce calibration plots", action='store_true')
    parser.add_argument("--doc_latex", help = "Generate a pdf document with latex", action='store_true')
    parser.add_argument("--export_odv4", help = "Export data in odv4 format", action='store_true')
    
    # parse the arguments
    args = parser.parse_args()

    smru_name = args.smru_name
    deployment = args.deployment
    
    
    if (smru_name or deployment) and \
        (args.do_all or \
         args.update_config_files or \
         args.import_data or \
         args.process_data or \
         args.create_hr2 or \
         args.calibration_plots or \
         args.doc_latex or \
         args.export_odv4):
        
        start_matlab()
        conf = init_mirounga()
        if args.update_config_files or args.do_all:
            conf = update_config_files(conf)
        if args.import_data or args.do_all:
            import_raw_data(deployment=deployment)
        if args.process_data or args.do_all:
            process_tags(deployment=deployment,smru_name=smru_name)
        if args.create_hr2 or args.do_all:
            create_hr2(deployment=deployment,smru_name=smru_name)
        if args.calibration_plots or args.do_all:
            generate_calibration_plots(deployment=deployment,smru_name=smru_name)
        if args.doc_latex or args.do_all:
            generate_doc_latex(deployment=deployment,smru_name=smru_name)
        if args.export_odv4 or args.do_all:
            export_odv4(deployment=deployment,smru_name=smru_name)
        stop_matlab()
        


            

