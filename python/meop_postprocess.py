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
import meop_filenames
import meop_metadata

processdir = meop_filenames.processdir


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


# flag profiles where LATITUDE=0
def remove_bad_locations():
    
    def flag_qc(fid,index,list_var):
        for var in list_var:
            if var in fid.variables:
                arr = fid[var][index]
                arr = np.where(arr==b'1',b'4',arr)
                fid[var][index] = arr
        return fid

    list_qf = ['lr0','lr1','hr0','hr1','fr0','fr1']
    list_var = ['TEMP_QC','TEMP_ADJUSTED_QC','PSAL_QC','PSAL_ADJUSTED_QC','CHLA_QC','CHLA_ADJUSTED_QC']
    
    for qf in list_qf:
        
        lprofiles, ltags, ldeployments = meop_metadata.read_lists_metadata(qf=qf)
        index = lprofiles[(lprofiles.LATITUDE==0)&(lprofiles.N_TEMP)].index
        if len(index)==0:
            continue
                
        list_smru_name = []
        for smru_name in lprofiles.loc[index,'SMRU_PLATFORM_CODE'].unique():
            namefile = meop_filenames.fname_prof(smru_name,qf=qf)
            if namefile.is_file():
                with netCDF4.Dataset(namefile,'a') as fid:
                    index2 = fid['LATITUDE'][:]==0
                    if index2.any():
                        flag_qc(fid,index2,list_var)
                        list_smru_name.append(smru_name)

        if list_smru_name:            
            file_pkl = processdir / f'list_profiles_{qf}.pkl'
            lprofiles = meop_metadata.update_lprofiles(lprofiles,list_smru_name,qf=qf)
            lprofiles.to_pickle(file_pkl)


def generate_descriptive_plots(smru_name='',deployment=''):
    
    list_qf = ['lr0','hr1','fr1']
    
    for qf in list_qf:

        for smru_name in meop.list_smru_name(smru_name,deployment,qf=qf):
            namefile = meop.fname_prof(smru_name,qf=qf)
            ds = meop.open_dataset(namefile)
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
    parser.add_argument("--metadata", help = "Update metadata", action='store_true')
    parser.add_argument("--descriptive_plots", help = "Produce descriptive plots", action='store_true')
    
    # parse the arguments
    args = parser.parse_args()

    smru_name = args.smru_name
    deployment = args.deployment
    
    
    if args.metadata or args.do_all:
        update_metadata(table_meta = table_meta, deployment=deployment,smru_name=smru_name)
    
    if args.descriptive_plots or args.do_all:
        generate_descriptive_plots(smru_name=smru_name,deployment=deployment)
        

        
