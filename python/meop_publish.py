from pathlib import Path
import os
import shutil
import xarray as xr
import pandas as pd
import numpy as np
import csv
import gsw
import matplotlib.pyplot as plt
from importlib import reload
import netCDF4 as nc

import meop
import meop_plot_data
import meop_metadata

# list functions

#  get_folder_public_name(config='jupyter_froqu_pc207_linux')
#  create_ncfile_all(smru_name,folder_out)
#  create_tag_plots(fname,folder_out,prefix_name)
#  publish_meop_ctd(folder_public, publish=True, genplots=True)

processdir = Path.home() / 'MEOP_process'

# copy a file
def copy_file(file_name,src_dir,dst_dir):
    shutil.copyfile(Path(src_dir)/file_name,Path(dst_dir)/file_name)

    
# get the name of folder where to publish in the congs.json file
def get_folder_public_name(config='jupyter_froqu_pc207_linux'):
    path = pd.read_json(processdir / 'configs.json').configs[config]['public']
    folder = pd.read_json(processdir / 'configs.json').version.CTDnew
    return Path(path) / folder
    

#  copy the variable var from nc_in in nc_out
def copy_netcdf_variable(nc_in,var_name_in,var_dims_in,nc_out,var_name_out,var_dims_out):
    
    with nc.Dataset(nc_in) as src, nc.Dataset(nc_out, "a") as dst:
        # copy dimensions if not already existing
        for i, name in enumerate(var_dims_out):
            if name not in dst.dimensions:
                dst.createDimension( name, src.dimensions[var_dims_in[i]].size )
            if src.dimensions[var_dims_in[i]].size - dst.dimensions[var_dims_out[i]].size != 0:
                print(f"Dimension {name} has wrong size in {nc_out}")
                return 0
        # copy variable
        if var_name_out not in dst.variables:
            var = dst.createVariable(var_name_out, src.variables[var_name_in].datatype, var_dims_out)
        # copy variable attributes all at once via dictionary
        dst[var_name_out].setncatts(src[var_name_in].__dict__)
        dst[var_name_out][:] = src[var_name_in][:]

    return 1


# create a netcdf file combining data on original levels and interp data
def create_ncfile_all(smru_name,folder_out):
    
    # copy lr1 in all
    ncfile_in = meop.fname_prof(smru_name,qf='lr1')
    ncfile_out = folder_out / meop.fname_prof(smru_name,qf='all').name

    if ncfile_in.is_file():
        
        shutil.copyfile(ncfile_in,ncfile_out)

        # copy ADJUSTED values in hr1 in INTERP variables
        ncfile_add = meop.fname_prof(smru_name,qf='hr1')
        copy_netcdf_variable(ncfile_add,'PRES_ADJUSTED',('N_PROF', 'N_LEVELS'),\
                                  ncfile_out,'PRES_INTERP',('N_PROF', 'N_INTERP'))
        copy_netcdf_variable(ncfile_add,'TEMP_ADJUSTED',('N_PROF', 'N_LEVELS'),\
                                  ncfile_out,'TEMP_INTERP',('N_PROF', 'N_INTERP'))
        copy_netcdf_variable(ncfile_add,'PSAL_ADJUSTED',('N_PROF', 'N_LEVELS'),\
                                  ncfile_out,'PSAL_INTERP',('N_PROF', 'N_INTERP'))
    
    return ncfile_out


# publish meop-ctd data in 
def create_tag_plots(fname,folder_out,prefix_name,var_suffix):

    with meop.open_dataset(fname) as ds:
        namefig = folder_out / (prefix_name+'_data_description.png')
        if not namefig.exists():
            fig, ax = ds.plot_data_tags('_ADJUSTED',namefig=namefig)
            plt.close(fig=fig)

    return


# publish meop-ctd data in 
def copy_license(folder_out):

    shutil.copyfile((processdir / 'README_licenseODbl.txt'),(Path(folder_out) / 'README_licenseODbl.txt'))

    return


# publish meop-ctd data in 
def publish_meop_ctd(folder_public, publish=True, genplots=True):

    folder_public = Path(folder_public)
    folder_public.mkdir(parents=True, exist_ok=True)
    if len(os.listdir(folder_public)):
        print(f'Warning: the public directory where to store public data {folder_public} is not empty. Risk of data corruption.')
    
    # copy license information
    copy_license(folder_public)

    list_profiles, list_tags, list_deployments = meop_metadata.read_list_profiles(rebuild=False,verbose=False,public=True,Tdata=False)

    for COUNTRY in list_deployments.COUNTRY.unique():
        
        print(COUNTRY)
        folder_country = folder_public / COUNTRY
        folder_country.mkdir(parents=True, exist_ok=True)
        folder_data = folder_country / 'DATA'
        folder_data.mkdir(parents=True, exist_ok=True)
        folder_plots = folder_public / COUNTRY / 'PLOTS'
        folder_plots.mkdir(parents=True, exist_ok=True)
        
        list_profiles_country, list_tags_country, list_deployments_country = \
            meop_metadata.filter_country(COUNTRY, list_profiles, list_tags, list_deployments)
        
        if publish:
            
            for smru_name in list_tags_country.SMRU_PLATFORM_CODE.unique():
            
                # copy ncfile: 'fr1 'if exists, otherwise 'all' with both low res and interp data
                is_done = (folder_data / meop.fname_prof(smru_name,qf='fr1').name).is_file() or \
                    (folder_data / meop.fname_prof(smru_name,qf='all').name).is_file()
                if not is_done:
                    print('Publish:',tag)
                    fname = meop.fname_prof(smru_name,qf='fr1')
                    if fname.exists():
                        shutil.copyfile(fname,folder_data / fname.name)
                    else:
                        fname = create_ncfile_all(tag,folder_data)
            
        if genplots:

            for smru_name in list_tags_country.SMRU_PLATFORM_CODE.unique():
            
                # create a plot for the tag
                namefig = folder_plots / (smru_name+'_data_description.png')
                if not namefig.is_file():
                    print('Generate plot:',smru_name)
                    # figure based on fr1 if possible. Otherwise based on adjusted profiles
                    fname = folder_data / meop.fname_prof(smru_name,qf='fr1').name
                    if fname.is_file():
                        create_tag_plots(fname,folder_plots,smru_name,'_ADJUSTED')
                    else:
                        fname = folder_data / meop.fname_prof(smru_name,qf='all').name
                        if fname.is_file():
                            create_tag_plots(fname,folder_plots,smru_name,'_INTERP')
        
    return



# Execute in terminal command line
if __name__ == "__main__":

    import argparse
    
    # create parser
    parser = argparse.ArgumentParser()

    # add arguments to the parser
    parser.add_argument("--path_public", help = "Provide path to public folder")
    parser.add_argument("--config_name", help = "Provide config name (see configs.json for a list). Not used if PATH_PUBLIC is provided.")
    parser.add_argument("--publish", help = "Publish data", action='store_true')
    parser.add_argument("--genplots", help = "Generate plots", action='store_true')
    
    # parse the arguments
    args = parser.parse_args()

    if args.path_public:
        folder_public = args.path_public
    else:
        if args.config_name:
            folder_public = get_folder_public_name(args.config_name)
        else:
            folder_public = get_folder_public_name()
            
    publish_meop_ctd(folder_public, publish=args.publish, genplots=args.genplots)
    
