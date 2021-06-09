from pathlib import Path
import os
import shutil
import xarray as xr
import pandas as pd
import numpy as np
import meop
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import gsw
from pathlib import Path
import cmocean.cm as cmo

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


# list variables
#  processdir
#
# list functions
#  deployment_from_smru_name(smru_name)
#  fname_prof(smru_name,deployment='',qf='lr0')
#  list_fname_prof(smru_name='',deployment='',qf='*')
#  list_smru_name(smru_name='',deployment='',qf='*')
#  fname_plots(smru_name,deployment='',qf='lr0',suffix='_plot')
#  N_PARAM(ds,PARAM)
#  copy_file(file_name,src_dir,dst_dir)
#  read_ncfile(ncfile_name)
#
# list of methods added to xr.Dataset objets loaded from ncARGO file
#  Decorator: @add_method(xr.Dataset)
#  central_longitude(self)
#  plot_map(self,ax=None,namefig=None,title='',figsize=(10, 10))
#  plot_profiles(self,PARAM,SUFFIX_PARAM='_ADJUSTED',ax=None,namefig=None,figsize=(10, 10))
#  plot_TSdiag(self,SUFFIX_PARAM='_ADJUSTED',ax=None,namefig=None,figsize=(10, 10))
#  plot_section(self,PARAM,SUFFIX_PARAM='_ADJUSTED',ax=None,namefig=None,figsize=(10, 10),title=None,rolling=0,**kwargs)
#  plot_data_tags(self,SUFFIX_PARAM='_ADJUSTED',namefig=None)
#  plot_TSsections(self,SUFFIX_PARAM='_ADJUSTED',namefig=None,figsize=(12, 10),rolling=0,**kwargs)



processdir = Path.home() / 'MEOP_process'

#-----------------------------------   utils     --------------------------------------------#

# 1. utils to reconstruct name of ncfiles
def deployment_from_smru_name(smru_name):
    return smru_name.split("-")[0]

# return ncARGO filename corresponding to a smru_name
def fname_prof(smru_name,deployment='',qf='lr0'):
    if not deployment:
        deployment = deployment_from_smru_name(smru_name)
    return Path(processdir,'final_dataset_prof',deployment,smru_name+'_'+qf+'_prof.nc')

# return list of ncARGO filename
def list_fname_prof(smru_name='',deployment='',qf='*'):
    if not deployment:
        deployment = deployment_from_smru_name(smru_name)
    dirEXP = Path(processdir,'final_dataset_prof',deployment)
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




# 2. utils to read list of deployments
# read list_deployment.csv file in processdir and return pandas dataframe
def read_list_deployment(filename=(processdir/'list_deployment.csv')):    
    if filename.is_file():
        list_deployment = pd.read_csv(filename)
        newnames = {}
        for var in list_deployment:
            newnames[var] = var.upper()
        return list_deployment.rename(columns=newnames).set_index('DEPLOYMENT_CODE')
    print('File',filename,'not found')
    return []

# read list_deployment.csv file in processdir and return pandas dataframe
def read_list_deployment_hr(filename=(processdir / 'list_deployment_hr.csv')):
    if filename.is_file():
        return pd.read_csv(filename, dtype={'prefix': str,'instr_id':str,'year':str})
    print('File',filename,'not found')
    return []

# return a DaraArray with the number of valid profile for the given PARAM
def N_PARAM(ds,PARAM):
    if PARAM+'_QC' in list(ds.variables):
        N_PARAM = np.sum(ds[PARAM+'_QC'].isin([b'1',b'8']),axis=1)
        N_PARAM = N_PARAM.where(N_PARAM!=0,np.nan)
    else:
        N_PARAM = xr.DataArray(np.empty(ds.dims['N_PROF']), dims=['N_PROF'])
        N_PARAM = np.nan
    return N_PARAM

# copy a file
def copy_file(file_name,src_dir,dst_dir):
    shutil.copyfile(Path(src_dir)/file_name,Path(dst_dir)/file_name)


    
    

    
    
#-----------------------------------  read ncARGO file  -------------------------------------#    
# read a netCDF file and return a xarray dataset structure
def read_ncfile(ncfile_name):
    
    ncfile_name = Path(ncfile_name)
    if ncfile_name.is_file():
        ds = xr.open_dataset(ncfile_name)
        for dim in ds.dims:
            ds[dim] = ((dim), ds[dim])
            ds.set_coords([dim])
    else:
        print('No file: ',ncfile_name)
        return None
    return ds
    


    

#----------------  add methods to xr.Dataset object with ncARGO format  -----------------------#    
# decorator to add a method to an object   
from functools import wraps
def add_method(cls):
    def decorator(func):
        @wraps(func) 
        def wrapper(self, *args, **kwargs): 
            return func(self, *args, **kwargs)
        setattr(cls, func.__name__, wrapper)
        # Note we are not binding func, but wrapper which accepts self but does exactly the same as func
        return func # returning func means func can still be used normally
    return decorator    
    
    

# compute number of valid profile by sensor type
@add_method(xr.Dataset)
def add_N_PARAM(self):
    ds = self
    ds['N_TEMP'] = N_PARAM(ds,'TEMP')
    ds['N_PSAL'] = N_PARAM(ds,'PSAL')
    if 'N_CHLA' in ds.variables:
        ds['N_CHLA'] = N_PARAM(ds,'CHLA')
    if 'N_DOXY' in ds.variables:
        ds['N_DOXY'] = N_PARAM(ds,'DOXY')
    return ds


            
# compute density and append to the dataset
@add_method(xr.Dataset)
def add_sigma0(self,SUFFIX_PARAM='_ADJUSTED'):    
    ds = self    
    if ('SIG0'+SUFFIX_PARAM) not in ds.variables:
        ds['SIG0'+SUFFIX_PARAM] = gsw.sigma0(ds['PSAL'+SUFFIX_PARAM],ds['TEMP'+SUFFIX_PARAM])
    return ds


# build list metadata in MEOP and return as dataframe
@add_method(xr.Dataset)
def list_metadata(self):

    ds=self
    data = {
        'DEPLOYMENT_CODE': deployment_from_smru_name(ds.smru_platform_code),
        'SMRU_PLATFORM_CODE': ds.smru_platform_code,
        'CYCLE_NUMBER': ds['CYCLE_NUMBER'].astype(int),
        'JULD': ds['JULD'],
        'LATITUDE': ds['LATITUDE'],
        'LONGITUDE': ds['LONGITUDE'],
        'N_TEMP' : meop.N_PARAM(ds,'TEMP'),
        'N_PSAL' : meop.N_PARAM(ds,'PSAL'),
        'N_CHLA' : meop.N_PARAM(ds,'CHLA')}
    df = pd.DataFrame(data)
    return df
            

# compute mixed layer depth and append to the dataframe
@add_method(xr.Dataset)
def add_mld(self,SUFFIX_PARAM='_ADJUSTED',density_threshold=0.02):    

    ds=self
    ds = ds.add_sigma0(SUFFIX_PARAM=SUFFIX_PARAM)
    density = ds['SIG0'+SUFFIX_PARAM].bfill(dim='N_LEVELS',limit=50)
    dens10 = density[:,9]
    pressure = ds.PRES.where(density-dens10<density_threshold,np.nan)
    mld = pressure.max(dim='N_LEVELS')
    mld[mld<5] = np.nan
    ds['MLD'+SUFFIX_PARAM] = mld
    return ds
            
            
            
            

            
            
            
# compute the value of central_longitude in order to create a map, based on the distribution of longitudes
@add_method(xr.Dataset)
def central_longitude(self):
    ds=self
    lon = ds['LONGITUDE']        
    if (lon.max()-lon.min()>180) & (np.abs(lon)>90).all():
        return 180
    else:
        return 0
        

        
# plot map
@add_method(xr.Dataset)
def plot_map(self,ax=None,namefig=None,title='',figsize=(10, 10),draw_background=True,color=[],scatter_plot=True):

    ds = self

    if not ax:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1,1,1,projection=ccrs.PlateCarree(central_longitude=ds.central_longitude()))
    else:
        fig = ax.get_figure()
    
    if scatter_plot:
        if list(color):
            h = ds.plot.scatter(ax=ax,x='LONGITUDE',y='LATITUDE',
                          color=color,add_guide=False,transform=ccrs.PlateCarree())
        else:
            h = ds.plot.scatter(ax=ax,x='LONGITUDE',y='LATITUDE',hue='N_PROF',hue_style='continuous',
                          add_guide=False,transform=ccrs.PlateCarree())
    else:
        if list(color):
            ax.plot(ds['LONGITUDE'].T,ds['LATITUDE'].T,linewidth=.6, color=color)
        else:
            ax.plot(ds['LONGITUDE'].T,ds['LATITUDE'].T,linewidth=.6)
        
        
    if draw_background:
        gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
          linewidth=1, color='gray', alpha=0.5)
        gl.top_labels = False
        gl.right_labels = False

        extent = list(ax.get_extent())
        if extent[1]-extent[0]<2:
            extent[0] -= 1.
            extent[1] += 1.
        if extent[3]-extent[2]<2:
            extent[2] -= 1.
            extent[3] += 1.
        ax.set_extent(extent, crs=ccrs.PlateCarree())
        #ax.stock_img()
        ax.coastlines()
        ax.add_feature(cfeature.LAND, color="lightgrey", alpha=0.5)

    if title:
        ax.set_title(title)

    if namefig:
        plt.savefig(namefig,dpi=300,bbox_inches='tight')

    return fig, ax

# plot_profiles(self,PARAM,SUFFIX_PARAM='_ADJUSTED',ax=None,namefig=None,figsize=(10, 10),color=[])
# PARAM can be a string 'TEMP','PSAL','CHLA'... or a list of strings
@add_method(xr.Dataset)
def plot_profiles(self,PARAM,SUFFIX_PARAM='_ADJUSTED',ax=None,namefig=None,figsize=(10, 10),color=[],pmax=1000):
    
    ds = self
    if isinstance(PARAM,str):
        if not ax:
            fig, ax = plt.subplots(1,1,figsize=figsize)
        else:
            fig = ax.get_figure()

        if PARAM=='SIG0':
            ds.add_sigma0()

        if PARAM=='SIG0':
            number_profiles = np.sum((N_PARAM(ds,'TEMP'+SUFFIX_PARAM)*N_PARAM(ds,'PSAL'+SUFFIX_PARAM)).data>0)
        else:
            number_profiles = np.sum(N_PARAM(ds,PARAM+SUFFIX_PARAM).data>0)

        if (PARAM+SUFFIX_PARAM) not in ds.variables.keys():
            return fig, ax
        
        # test if color is fixed
        if list(color):
            ax.plot(ds[PARAM+SUFFIX_PARAM].T,ds['PRES'+SUFFIX_PARAM].T,linewidth=.6, color=color)
        else:
            ax.plot(ds[PARAM+SUFFIX_PARAM].T,ds['PRES'+SUFFIX_PARAM].T,linewidth=.6)
            
        if 'MLD_ADJUSTED' in ds.variables:
            ax.scatter(ds[PARAM+SUFFIX_PARAM].sel(N_LEVELS=ds.MLD_ADJUSTED),ds.MLD_ADJUSTED,color='k')
            
        ax.set_ylim([pmax,0])

        ax.set_title(f"{PARAM+SUFFIX_PARAM}: {number_profiles} profiles")

        if namefig:
            plt.savefig(namefig,dpi=300,bbox_inches='tight')

        return fig, ax
    
    elif isinstance(PARAM,list):
        
        if len(PARAM)==1:
            ds.plot_profiles(PARAM[0],SUFFIX_PARAM=SUFFIX_PARAM,ax=ax,namefig=namefig,figsize=figsize,color=color,pmax=pmax)
        else:
            if not ax:
                fig, ax = plt.subplots(1,len(PARAM),figsize=figsize)
            else:
                fig = ax.get_figure()
            for kk,param in enumerate(PARAM):
                ds.plot_profiles(param,SUFFIX_PARAM=SUFFIX_PARAM,ax=ax[kk],color=color,pmax=pmax)
            if namefig:
                plt.savefig(namefig,dpi=300,bbox_inches='tight')

            return fig, ax

    


@add_method(xr.Dataset)
def plot_TSdiag(self,SUFFIX_PARAM='_ADJUSTED',mode='line',ax=None,namefig=None,figsize=(10, 10),draw_sigma=True,color=[]):

    ds = self

    if not ax:
        fig, ax = plt.subplots(1,1,figsize=figsize)
    else:
        fig = ax.get_figure()
        
    # add pressure
    if SUFFIX_PARAM == '_INTERP':
        name_coord_levels = "N_INTERP"
        ds = ds.assign_coords(pressure=(name_coord_levels, ds.PRES_INTERP[0,:]))
    else: # processing = {'','_ADJUSTED'}
        name_coord_levels = "N_LEVELS"
        ds = ds.assign_coords(pressure=(name_coord_levels, ds.PRES[0,:]))
    
    if list(color):
        ax.plot(ds['PSAL'+SUFFIX_PARAM].T,ds['TEMP'+SUFFIX_PARAM].T,linewidth=.6,color=color)
    else:
        ax.plot(ds['PSAL'+SUFFIX_PARAM].T,ds['TEMP'+SUFFIX_PARAM].T,linewidth=.6)
                
    if draw_sigma:
        x0,x1 = ax.get_xlim()
        y0,y1 = ax.get_ylim()
        X,Y = np.meshgrid(np.linspace(x0,x1),np.linspace(y0,y1))
        D = gsw.sigma0(X,Y)
        ax.set_aspect((x1-x0)/(y1-y0))
        ax.contour(X,Y,D,colors='k')
    
    if namefig:
        plt.savefig(namefig,dpi=300,bbox_inches='tight')

    return fig, ax


@add_method(xr.Dataset)
def plot_sections(self,PARAM,SUFFIX_PARAM='_ADJUSTED',ax=None,namefig=None,figsize=(10, 10),title=None,rolling=0,pmax=1000,**kwargs):


    ds = self

    if isinstance(PARAM,str):
        
        if not ax:
            fig, ax = plt.subplots(1,1,figsize=figsize)
        else:
            fig = ax.get_figure()

        if (PARAM+SUFFIX_PARAM) not in ds.variables.keys():
            return fig, ax

        if 'cmap' not in kwargs:
            if PARAM=='TEMP':
                kwargs['cmap'] = cmo.thermal
            elif PARAM=='PSAL':
                kwargs['cmap'] = cmo.haline
            else:
                kwargs['cmap'] = cmo.solar
        if 'yincrease' not in kwargs:
            kwargs['yincrease'] = False

        if rolling:
            var = ds[PARAM+SUFFIX_PARAM].rolling(N_PROF=rolling, center=True, min_periods=1).mean()
        else:
            var = ds[PARAM+SUFFIX_PARAM]
        var.T.plot(ax=ax,**kwargs)    

        if not title:
            number_profiles = np.sum(N_PARAM(ds,PARAM+SUFFIX_PARAM).data>0)
            title = f"{ds.smru_platform_code}, {PARAM+SUFFIX_PARAM}: {number_profiles} profiles"
        ax.set_ylim(pmax)
        ax.set_title(title)
        
        # save figure
        if namefig:
            plt.savefig(namefig,dpi=300)

        return fig, ax

    elif isinstance(PARAM,list):
        
        if len(PARAM)==1:
            ds.plot_sections(PARAM[0],SUFFIX_PARAM=SUFFIX_PARAM,ax=ax,namefig=namefig,figsize=figsize,title=title,rolling=rolling,pmax=pmax,**kwargs)
        else:
            if not ax:
                fig, ax = plt.subplots(len(PARAM),1,figsize=figsize)
            else:
                fig = ax.get_figure()
            for kk,param in enumerate(PARAM):
                ds.plot_sections(param,SUFFIX_PARAM=SUFFIX_PARAM,ax=ax[kk],title=title,rolling=rolling,pmax=pmax,**kwargs)
            if namefig:
                plt.savefig(namefig,dpi=300,bbox_inches='tight')

            return fig, ax

        
@add_method(xr.Dataset)
def plot_TSsections(self,SUFFIX_PARAM='_ADJUSTED',ax=None,namefig=None,figsize=(10, 10),title=None,rolling=0,pmax=1000,**kwargs):

    ds = self
    
    fig, ax = ds.plot_sections(['TEMP','PSAL'],\
                                 SUFFIX_PARAM=SUFFIX_PARAM,ax=ax,namefig=namefig,figsize=figsize,title=title,pmax=pmax,rolling=rolling,**kwargs)
    if SUFFIX_PARAM=='_ADJUSTED':
        ds.add_mld()
        mld = ds.MLD_ADJUSTED.rolling(N_PROF=rolling,min_periods=3,center=True).mean()
        for a in ax:
            mld.plot(ax=a)
    return fig, ax


        
        
@add_method(xr.Dataset)
def plot_data_tags(self,SUFFIX_PARAM='_ADJUSTED',namefig=None):

    ds = self
    
    if ds['TEMP'+SUFFIX_PARAM].sum(dim='N_LEVELS').sum().data==0:
        return None, None
    
    from cycler import cycler
    cmap = plt.get_cmap('viridis',ds.dims['N_PROF'])
    custom_cycler = cycler(color=cmap.colors)

    fig = plt.figure(figsize=(12,9))
    gs = fig.add_gridspec(2,3)

    ax={}
    ax['T'] = fig.add_subplot(gs[1, 0])
    ax['S'] = fig.add_subplot(gs[1, 1])
    ax['D'] = fig.add_subplot(gs[1, 2])
    ax['TS'] = fig.add_subplot(gs[0, 2:])
    ax['xy'] = fig.add_subplot(gs[0, :2],projection=ccrs.PlateCarree(central_longitude=ds.central_longitude()))
    for key in ax:
        ax[key].set_prop_cycle(custom_cycler)

    ds.plot_profiles('TEMP',SUFFIX_PARAM,ax=ax['T'])
    ds.plot_profiles('PSAL',SUFFIX_PARAM,ax=ax['S'])
    ds.plot_profiles('SIG0',SUFFIX_PARAM,ax=ax['D'])
    ds.plot_TSdiag(SUFFIX_PARAM,ax=ax['TS'])

    # main title
    start_date = ds.JULD.min()
    end_date = ds.JULD.max()
    num_days = (end_date-start_date).dt.days
    title = f"{ds.smru_platform_code}: from {start_date.dt.strftime('%Y-%m-%d').data} to {end_date.dt.strftime('%Y-%m-%d').data} ({num_days.data} days)"
    
    # map
    ds.plot_map(ax=ax['xy'], title=title)

    # finitions
    plt.tight_layout()

    # save figure
    if namefig:
        plt.savefig(namefig,dpi=300)

    return fig, ax



