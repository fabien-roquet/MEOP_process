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
import meop_filenames

# list functions

#  create_ncfile_all(smru_name,folder_out)
#  create_tag_plots(fname,folder_out,prefix_name)
#  publish_meop_ctd(folder_public, publish=True, genplots=True)

processdir = meop_filenames.processdir

# copy a file
def copy_file(file_name,src_dir,dst_dir):
    shutil.copyfile(Path(src_dir)/file_name,Path(dst_dir)/file_name)

    
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
    ncfile_in = meop_filenames.fname_prof(smru_name,qf='lr1')
    ncfile_out = folder_out / meop_filenames.fname_prof(smru_name,qf='all').name

    if ncfile_in.is_file():
        
        shutil.copyfile(ncfile_in,ncfile_out)

        # copy ADJUSTED values in hr1 in INTERP variables
        ncfile_add = meop_filenames.fname_prof(smru_name,qf='hr1')
        copy_netcdf_variable(ncfile_add,'PRES_ADJUSTED',('N_PROF', 'N_LEVELS'),\
                                  ncfile_out,'PRES_INTERP',('N_PROF', 'N_INTERP'))
        copy_netcdf_variable(ncfile_add,'TEMP_ADJUSTED',('N_PROF', 'N_LEVELS'),\
                                  ncfile_out,'TEMP_INTERP',('N_PROF', 'N_INTERP'))
        copy_netcdf_variable(ncfile_add,'PSAL_ADJUSTED',('N_PROF', 'N_LEVELS'),\
                                  ncfile_out,'PSAL_INTERP',('N_PROF', 'N_INTERP'))
    
    return ncfile_out


# publish meop-ctd data in 
def copy_license(folder_out):
    shutil.copyfile((processdir / 'README_licenseODbl.txt'),(Path(folder_out) / 'README_licenseODbl.txt'))
    return


# list_deployment obtained with read_list_deployment()
def load_list_profiles(publicdir_CTD, rebuild=False):
    
    if rebuild or not (publicdir_CTD / 'list_profiles.csv').exists():

        # list of deployment dataframes
        list_df=[]

        # get list of metadata for each tag and store in a list
        for ncfile in publicdir_CTD.glob('**/*_all_prof.nc'):

            if ncfile.exists():
                with meop.open_dataset(ncfile) as ds:

                    csv_list_file = ncfile.parent / f"{ncfile.stem.split('_')[0]}_list.csv"
                    if not csv_list_file.exists():
                        df = ds.list_metadata()
                        df.to_csv(csv_list_file ,index=False)
                    else:
                        df = pd.read_csv(csv_list_file)
                list_df.append(df)

        if not list_df:
            return []

        # concatenate list of dataframes into one dataframe
        lprofiles = pd.concat(list_df)
        lprofiles.to_csv(publicdir_CTD / 'list_profiles.csv',index=False)

        ltags, ldeployments = meop_metadata.list_tags_deployments_from_profiles(lprofiles)
        ltags.to_csv(publicdir_CTD / 'list_tags.csv',index=False)
        ldeployments.to_csv(publicdir_CTD / 'list_deployments.csv',index=False)
       
    else:
        
        lprofiles = pd.read_csv(publicdir_CTD / 'list_profiles.csv')
        ltags, ldeployments = meop_metadata.list_tags_deployments_from_profiles(lprofiles)

    return lprofiles, ltags, ldeployments


# add maps
def build_maps(publicdir_CTD, rebuild = False):
    
    lprofiles, ltags, ldeployments = load_list_profiles(publicdir_CTD, rebuild=False)
    lprofiles = lprofiles[lprofiles.N_TEMP>0]
    ltags, ldeployments = meop_metadata.list_tags_deployments_from_profiles(lprofiles)
    if 'COUNTRY' not in lprofiles:
        lprofiles = lprofiles.merge(ldeployments[['DEPLOYMENT_CODE','COUNTRY']],on='DEPLOYMENT_CODE',right_index=False,)

    meop_plot_data.plot_map_deployments(lprofiles,
                              groupby='MASK',
                              title=f'Distribution of profiles',
                              legend=True, show_plot=False,
                              namefig='Global_distribution_by_regions.png',
                              folder=publicdir_CTD,
                             )    

    meop_plot_data.plot_map_deployments(lprofiles,
                              groupby='DEPLOYMENT_CODE',
                              title=f'Distribution of profiles by deployment code',
                              legend=False, show_plot=False,
                              namefig='Global_distribution_by_deployment.png',
                              folder=publicdir_CTD,
                             )

    meop_plot_data.plot_map_deployments(lprofiles,
                              groupby='DEPLOYMENT_CODE',
                              title=f'Distribution of profiles by deployment code',
                              legend=False, show_plot=False,
                              legend_horiz=True,
                              namefig='Global_distribution_by_deployment_with_legend.png',
                              folder=publicdir_CTD,
                             )    

    meop_plot_data.plot_map_deployments(lprofiles,
                              groupby='COUNTRY',
                              title=f'Distribution of profiles by country',
                              legend=True, show_plot=False,
                              folder=publicdir_CTD, namefig='Global_distribution_by_country.png',                              
                             )
    
    for region in list(ltags.MASK.unique()):
        index_tags = ltags[ltags.MASK==region].SMRU_PLATFORM_CODE
        meop_plot_data.plot_map_deployments(lprofiles[lprofiles.SMRU_PLATFORM_CODE.isin(index_tags)],
                                  groupby='DEPLOYMENT_CODE',
                                  title=f'Distribution of {region} profiles',
                                  legend=True, legend_horiz=True, show_plot=False,
                                  folder=publicdir_CTD, namefig=f'Regional_distribution_{region}.png',                                  
                                 )
    
    for country in list(lprofiles.COUNTRY.unique()):
        lprof_country = lprofiles[lprofiles.COUNTRY==country]
        meop_plot_data.plot_map_deployments(lprof_country,
                                  groupby='DEPLOYMENT_CODE',
                                  title=f'Distribution of profiles for {country}',
                                  legend=True, legend_horiz=True, show_plot=False,
                                  folder=publicdir_CTD, namefig=f'National_distribution_{country}.png',                                  
                                 )    

    for country in list(lprofiles.COUNTRY.unique()):
        ldepl_country = ldeployments[ldeployments.COUNTRY==country]
        folder = publicdir_CTD / country
        for depl in ldepl_country.DEPLOYMENT_CODE:
            # map for each deployment
            namefig= f"map_{depl}.png"
            if not (folder/namefig).exists() or rebuild:
                meop_plot_data.plot_map_deployments(lprofiles[lprofiles.DEPLOYMENT_CODE==depl],\
                                                   folder = folder, namefig=namefig, show_plot=False)
                
            # info for each deployment
            namefig= folder / f"info_{depl}.png"
            if not (folder/namefig).exists() or rebuild:
                list_ncfile = list(folder.glob(f'**/{depl}*_all_prof.nc'))
                meop_plot_data.plot_data_deployments(depl,namefig=namefig,list_fname_prof=list_ncfile)
                
    return


# publish meop-ctd data in 
def create_tag_plots(fname,folder_out,prefix_name,var_suffix):

    with meop.open_dataset(fname) as ds:
        namefig = folder_out / (prefix_name+'_data_description.png')
        if not namefig.exists():
            fig, ax = ds.plot_data_tags('_ADJUSTED',namefig=namefig)
            plt.close(fig=fig)

    return




# publish meop-ctd data in 
def publish_meop_ctd(folder_public, publish=True, genplots=True, genmaps=True, verbose=False):

    folder_public = Path(folder_public)
    folder_public.mkdir(parents=True, exist_ok=True)
    if len(os.listdir(folder_public)):
        print(f'Warning: the public directory where to store public data {folder_public} is not empty. Risk of data corruption.')
    
    # copy license information
    copy_license(folder_public)

    lprofiles, ltags, ldeployments = meop_metadata.read_lists_metadata(rebuild=False,verbose=False,public=True,Tdata=False)

    for COUNTRY in ldeployments.COUNTRY.unique():
        
        print(COUNTRY)
        folder_country = folder_public / COUNTRY
        folder_country.mkdir(parents=True, exist_ok=True)
        folder_data = folder_country / 'DATA'
        folder_data.mkdir(parents=True, exist_ok=True)
        folder_plots = folder_public / COUNTRY / 'PLOTS'
        folder_plots.mkdir(parents=True, exist_ok=True)
        
        lprofiles_country, ltags_country, ldeployments_country = \
            meop_metadata.filter_country(COUNTRY, lprofiles, ltags, ldeployments)
        
        if publish:
            
            print('...Generation of netCDF files')
            for smru_name in ltags_country.SMRU_PLATFORM_CODE.unique():
            
                # copy ncfile: 'fr1 'if exists, otherwise 'all' with both low res and interp data
                is_done = (folder_data / meop_filenames.fname_prof(smru_name,qf='fr1').name).is_file() or \
                    (folder_data / meop_filenames.fname_prof(smru_name,qf='all').name).is_file()
                if not is_done:
                    if verbose:
                        print('Publish:',smru_name)
                    fname = meop_filenames.fname_prof(smru_name,qf='fr1')
                    if fname.exists():
                        shutil.copyfile(fname,folder_data / fname.name)
                    else:
                        fname = create_ncfile_all(smru_name,folder_data)
            
        if genplots:

            print('...Generation of descriptive plots')
            for smru_name in ltags_country.SMRU_PLATFORM_CODE.unique():
            
                # create a plot for the tag
                namefig = folder_plots / (smru_name+'_data_description.png')                
                if not namefig.is_file():
                    if verbose:
                        print('Generate plot:',smru_name)
                    # figure based on fr1 if possible. Otherwise based on adjusted profiles
                    fname = folder_data / meop_filenames.fname_prof(smru_name,qf='fr1').name
                    if fname.is_file():
                        create_tag_plots(fname,folder_plots,smru_name,'_ADJUSTED')
                    else:
                        fname = folder_data / meop_filenames.fname_prof(smru_name,qf='all').name
                        if fname.is_file():
                            create_tag_plots(fname,folder_plots,smru_name,'_INTERP')
        
    if genmaps:

        print('Generation of maps')
        build_maps(folder_public)
            
    return



# Execute in terminal command line
if __name__ == "__main__":

    import argparse
    
    # create parser
    parser = argparse.ArgumentParser()

    # add arguments to the parser
    parser.add_argument("--path_public", help = "Provide path to public folder")
    parser.add_argument("--publish", help = "Publish data", action='store_true')
    parser.add_argument("--genplots", help = "Generate plots", action='store_true')
    parser.add_argument("--genmaps", help = "Generate maps", action='store_true')
    
    # parse the arguments
    args = parser.parse_args()

    if args.path_public:
        folder_public = Path(args.path_public)
    else:
        folder_public = meop_filenames.publicdir_CTD
    print('Publish in public folder: '+str(folder_public))
    
    publish_meop_ctd(folder_public, publish=args.publish, genplots=args.genplots, genmaps=args.genmaps)
    
