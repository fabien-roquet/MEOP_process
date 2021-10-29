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

    if ncfile_in.is_file() and not ncfile_out.is_file():
        
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
def load_list_profiles(publicdir_CTD, public=True, rebuild=False):
    
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
        
        # save lists
        lprofiles.to_csv(publicdir_CTD / 'list_profiles.csv',index=False)
        ltags, ldeployments = meop_metadata.list_tags_deployments_from_profiles(lprofiles)
        ltags.to_csv(publicdir_CTD / 'list_tags.csv',index=False)
        ldeployments.to_csv(publicdir_CTD / 'list_deployments.csv',index=False)

    else:
        
        lprofiles = pd.read_csv(publicdir_CTD / 'list_profiles.csv')
        
    if public:
        lprofiles = lprofiles[lprofiles.N_TEMP>0]
        
    ltags, ldeployments = meop_metadata.list_tags_deployments_from_profiles(lprofiles)
       
    if 'COUNTRY' not in lprofiles:
        lprofiles = lprofiles.merge(ldeployments[['DEPLOYMENT_CODE','COUNTRY']],on='DEPLOYMENT_CODE',right_index=False,)

    return lprofiles, ltags, ldeployments


# add plots
def build_plots(publicdir_CTD, rebuild = False, verbose=False):

    lprofiles, ltags, ldeployments = load_list_profiles(publicdir_CTD, public=True, rebuild=False)
    
    for country in list(lprofiles.COUNTRY.unique()):
        
        ldepl_country = ldeployments[ldeployments.COUNTRY==country]
        folder_country = publicdir_CTD / country
        folder_data = folder_country / 'DATA'
        folder_plots = folder_country / 'PLOTS'
        
        for depl in ldepl_country.DEPLOYMENT_CODE:

            # info for each deployment
            namefig= folder_country / f"info_{depl}.png"
            if not (folder_country/namefig).exists() or rebuild:
                list_ncfile = list(folder_country.glob(f'**/{depl}*_all_prof.nc'))
                meop_plot_data.plot_data_deployments(depl,namefig=namefig,list_fname_prof=list_ncfile)

            # info for each tag
            lprof = lprofiles[lprofiles.DEPLOYMENT_CODE==depl]
            for smru_name in list(lprof.SMRU_PLATFORM_CODE.unique()):

                # create a plot for the tag
                namefig = folder_plots / (smru_name+'_data_description.png')                
                if not namefig.is_file() or rebuild:
                    if verbose:
                        print('Generate plot:',smru_name)
                    # figure based on fr1 if possible. Otherwise based on adjusted profiles
                    fname = folder_data / meop_filenames.fname_prof(smru_name,qf='all').name
                    if fname.is_file():
                        create_tag_plots(fname,folder_plots,smru_name,'_INTERP')

            # info for each tag
            lprof = lprofiles[lprofiles.DEPLOYMENT_CODE==depl]
            for smru_name in list(lprof.SMRU_PLATFORM_CODE.unique()):

                # figure based on all-INTERP. Otherwise based on adjusted profiles
                fname = folder_data / meop_filenames.fname_prof(smru_name,qf='all').name
                if fname.is_file():
                    with meop.open_dataset(fname) as ds:
                        namefig = folder_plots / (smru_name+'_data_description.png')
                        if not namefig.is_file() or rebuild:
                            ds.plot_data_tags('_INTERP',namefig=namefig)
                    with meop.open_dataset(fname) as ds:
                        namefig = folder_plots / (smru_name+'_data_sections.png')   
                        if not namefig.is_file() or rebuild:
                            ds.plot_sections(['TEMP','PSAL','SIG0'],rolling=1,density_threshold=0.03,namefig=namefig)

    return


# add maps
def build_maps(publicdir_CTD, rebuild = False):
    
    lprofiles, ltags, ldeployments = load_list_profiles(publicdir_CTD, public=True, rebuild=False)

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
                                
    return



# compress files. Stored in the parent directory of publicdir_CTD
def compress_public_data(publicdir_CTD, rebuild = False):

    
    lprofiles, ltags, ldeployments = load_list_profiles(publicdir_CTD, public=True, rebuild=False)
    version = publicdir_CTD.name
    
    for COUNTRY in ldeployments.COUNTRY.unique():
        zip_file = publicdir_CTD.parent / f"{version}_{COUNTRY}.zip"
        if zip_file.exists() and rebuild:
            os.remove(zip_file)
        if (publicdir_CTD / COUNTRY).exists() and not zip_file.exists():
            bashCommand = f"zip -r {zip_file} {publicdir_CTD / COUNTRY}"
            print(bashCommand)
            os.system(bashCommand)
            
    zip_file = publicdir_CTD.parent / f"{version}_ALL.zip"
    if zip_file.exists() and rebuild:
        os.remove(zip_file)
    if not zip_file.exists():
        bashCommand = f"zip -r {zip_file} {publicdir_CTD}"
        print(bashCommand)
        os.system(bashCommand)

    return
    

# publish meop-ctd data in 
def publish_meop_ctd(publicdir_CTD=meop_filenames.publicdir_CTD, copydata=False, genplots=False, genmaps=False, compress=False, rebuild=False, verbose=False):

    publicdir_CTD.mkdir(parents=True, exist_ok=True)
    if len(os.listdir(publicdir_CTD)):
        print(f'Warning: the public directory where to store public data {publicdir_CTD} is not empty. Risk of data corruption.')
    
    # copy license information
    copy_license(publicdir_CTD)

    if copydata:

        print('...Generation of netCDF files')

        lprofiles, ltags, ldeployments = meop_metadata.read_lists_metadata(rebuild=False,verbose=False,public=True,Tdata=False)

        for COUNTRY in ldeployments.COUNTRY.unique():
            folder_country = publicdir_CTD / COUNTRY
            folder_country.mkdir(parents=True, exist_ok=True)
            folder_data = folder_country / 'DATA'
            folder_data.mkdir(parents=True, exist_ok=True)
            folder_plots = folder_country / 'PLOTS'
            folder_plots.mkdir(parents=True, exist_ok=True)
    
        for COUNTRY in ldeployments.COUNTRY.unique():

            print(COUNTRY)
            folder_country = publicdir_CTD / COUNTRY
            folder_data = folder_country / 'DATA'
            folder_plots = folder_country / 'PLOTS'

            lprofiles_country, ltags_country, ldeployments_country = \
                meop_metadata.filter_country(COUNTRY, lprofiles, ltags, ldeployments)

            for smru_name in ltags_country.SMRU_PLATFORM_CODE.unique():

                # copy ncfile: 'fr1' if exists
                folder_data_fr = folder_country / 'DATA_FULL_RES'
                fname_orig = meop_filenames.fname_prof(smru_name,qf='fr1')
                fname_copy = folder_data_fr / fname_orig.name
                if rebuild and fname_copy.is_file():
                    os.remove(fname_copy)
                if not fname_copy.is_file() and fname_orig.is_file():
                    if verbose:
                        print('Publish: ',fname_orig.name)          
                    folder_data_fr.mkdir(parents=True, exist_ok=True)                        
                    shutil.copyfile(fname_orig,fname_copy)
                    
                # copy ncfile: create 'all' from hr1 and lr1
                fname_copy = folder_data / meop_filenames.fname_prof(smru_name,qf='all').name
                if rebuild and fname_copy.is_file():
                    os.remove(fname_copy)
                if not fname_copy.is_file():
                    if verbose:
                        print('Publish: ',fname_orig.name)
                    create_ncfile_all(smru_name,folder_data)
                            
    if genplots:

        print('...Generation of plots')
        build_plots(publicdir_CTD,rebuild=rebuild)
        
        
    if genmaps:

        print('...Generation of maps')
        build_maps(publicdir_CTD,rebuild=rebuild)
            
    if compress:
        
        print('...Compress files')
        compress_public_data(publicdir_CTD,rebuild=rebuild)
        
    return



# Execute in terminal command line
if __name__ == "__main__":

    import argparse
    
    # create parser
    parser = argparse.ArgumentParser()

    # add arguments to the parser
    parser.add_argument("--path_public", help = "Provide path to public folder")
    parser.add_argument("--version", help = "Provide version of database")
    parser.add_argument("--copydata", help = "Copy data", action='store_true')
    parser.add_argument("--genplots", help = "Generate plots", action='store_true')
    parser.add_argument("--genmaps", help = "Generate maps", action='store_true')
    parser.add_argument("--compress", help = "Compress files", action='store_true')
    parser.add_argument("--rebuild", help = "delete and replace previous version of file", action='store_true')
    
    # parse the arguments
    args = parser.parse_args()

    if args.path_public:
        folder_public = Path(args.path_public)
    else:
        folder_public = meop_filenames.publicdir
        
    if args.version:
        version = args.version
    else:
        version = meop_filenames.version
 
    publicdir_CTD = folder_public / version
    print('Publish in public folder: '+str(publicdir_CTD))
    
    publish_meop_ctd(publicdir_CTD, copydata=args.copydata, genplots=args.genplots, genmaps=args.genmaps, compress=args.compress, rebuild=args.rebuild)
    
