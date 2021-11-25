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



    


# Execute in terminal command line
if __name__ == "__main__":

    import argparse
    
    # create parser
    parser = argparse.ArgumentParser()

    # add arguments to the parser
    parser.add_argument("--do_all", help = "Process data and produce plots", action='store_true')
    parser.add_argument("--remove_bad_locations", help = "Produce descriptive plots", action='store_true')
    
    # parse the arguments
    args = parser.parse_args()
    
    if args.remove_bad_locations or args.do_all:
        remove_bad_locations()
        

        
