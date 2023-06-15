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
import cmocean.cm as cmo

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)



    
#-----------------------------------  read ncARGO file  -------------------------------------#    
# read a netCDF file and return a xarray dataset structure
def open_dataset(ncfile_name):
    
    ncfile_name = Path(ncfile_name)
    if ncfile_name.is_file():
        ds = xr.open_dataset(ncfile_name)
        for dim in ds.dims:
            ds[dim] = xr.DataArray( data = ds[dim], dims = [dim])
            ds.set_coords(dim)
        ds['JULD'].values = [dt.replace(microsecond=0)  for dt in ds['JULD'].values]
        ds['JULD_LOCATION'] = ds.JULD
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
def N_PARAM(ds,PARAM,SUFFIX_PARAM='_ADJUSTED'):
    if PARAM+SUFFIX_PARAM+'_QC' in list(ds.variables):
        N_PARAM = np.sum(ds[PARAM+SUFFIX_PARAM+'_QC'].isin([b'1',b'8']),axis=1)
        N_PARAM = N_PARAM.where(N_PARAM!=0,np.nan)
    elif PARAM+SUFFIX_PARAM in list(ds.variables):
        N_PARAM = (~ds[PARAM+SUFFIX_PARAM].isnull()).sum(axis=1)
        N_PARAM = N_PARAM.where(N_PARAM!=0,np.nan)
    else:
        # filled with NaNs
        N_PARAM = xr.DataArray(dims=['N_PROF'],coords={'N_PROF':np.arange(ds.dims['N_PROF'])})
    return N_PARAM



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


            
# compute sigma0 and append to the dataset (SIG0 or SIG0_ADJUSTED), as well as its corresponding QC
@add_method(xr.Dataset)
def add_sigma0(self,SUFFIX_PARAM='_ADJUSTED'):    
    ds = self
    if ('SIG0'+SUFFIX_PARAM) not in ds.variables:
        if SUFFIX_PARAM=='_INTERP':
            ds = ds.add_interp('TEMP')
            ds = ds.add_interp('PSAL')
        ds['SIG0'+SUFFIX_PARAM] = gsw.sigma0(ds['PSAL'+SUFFIX_PARAM],ds['TEMP'+SUFFIX_PARAM])
        if not SUFFIX_PARAM=='_INTERP':
            ds['SIG0'+SUFFIX_PARAM+'_QC'] = ds['TEMP'+SUFFIX_PARAM+'_QC'].copy()
            ds['SIG0'+SUFFIX_PARAM+'_QC'].where(~ds['PSAL'+SUFFIX_PARAM+'_QC'].isin([b'4']),'4')
    return ds


# compute mixed layer depth and append to the dataframe
def compute_mld(ds,SUFFIX_PARAM='_ADJUSTED',density_threshold=0.02):
    N_LEVELS = 'N_LEVELS'
    if SUFFIX_PARAM=='_INTERP':
        N_LEVELS = 'N_INTERP'
    ds = ds.add_sigma0(SUFFIX_PARAM=SUFFIX_PARAM)
    density = ds['SIG0'+SUFFIX_PARAM].bfill(dim=N_LEVELS,limit=50)
    dens10 = density[:,9]
    pressure = ds['PRES'+SUFFIX_PARAM].where(density-dens10<density_threshold,np.nan)
    mld = pressure.max(dim=N_LEVELS)
    mld[mld<5] = np.nan
    return mld
            
            
# compute mixed layer depth and append to the dataframe
@add_method(xr.Dataset)
def add_mld(self,SUFFIX_PARAM='_ADJUSTED',density_threshold=0.02):    

    ds = self
    mld = compute_mld(ds,SUFFIX_PARAM=SUFFIX_PARAM,density_threshold=density_threshold)
    ds['MLD'+SUFFIX_PARAM] = mld
    return ds

            
            

# interpolate data pon a regular verticl grid and append to the dataframe
@add_method(xr.Dataset)
def add_interp(self,PARAM,SUFFIX_PARAM='_ADJUSTED',pmax=1000):
    ds = self
    if PARAM+'_INTERP' not in ds.variables:
        if (PARAM+'_ADJUSTED') not in ds.variables.keys():
            return ds
        N_INTERP = pmax
        data = np.zeros([ds.dims['N_PROF'],N_INTERP])*np.nan
        ds = ds.assign_coords(dict(N_INTERP=("N_INTERP", np.arange(0,N_INTERP))))
        da = xr.DataArray(data,dims=["N_PROF", "N_INTERP"],coords=dict(N_PROF=("N_PROF", ds.coords['N_PROF'].values),N_INTERP=("N_INTERP", ds.coords['N_INTERP'].values)))
        for i in range(ds.dims['N_PROF']):
            da[i,:] = np.interp(da.N_INTERP,ds['PRES_ADJUSTED'].values[i,:],ds[PARAM+'_ADJUSTED'].values[i,:])
        ds[PARAM+'_INTERP'] = da
        if 'PRES_INTERP' not in ds.variables:
            ds['PRES_INTERP'] = np.arange(0,N_INTERP)
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
            ax.scatter(ds.LONGITUDE,ds.LATITUDE,s=5,c=color,transform=ccrs.PlateCarree())
        else:
            ax.scatter(ds.LONGITUDE,ds.LATITUDE,s=5,c=ds.N_PROF,transform=ccrs.PlateCarree())
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
            ds.add_sigma0(SUFFIX_PARAM=SUFFIX_PARAM)

        if PARAM=='SIG0':
            number_profiles = np.sum((N_PARAM(ds,'TEMP',SUFFIX_PARAM)*N_PARAM(ds,'PSAL',SUFFIX_PARAM)).data>0)
        else:
            number_profiles = np.sum(N_PARAM(ds,PARAM,SUFFIX_PARAM).data>0)

        if (PARAM+SUFFIX_PARAM) not in ds.variables.keys():
            return fig, ax
        
        # test if color is fixed
        if list(color):
            ax.plot(ds[PARAM+SUFFIX_PARAM].T,ds['PRES'+SUFFIX_PARAM].T,linewidth=.6, color=color)
        else:
            ax.plot(ds[PARAM+SUFFIX_PARAM].T,ds['PRES'+SUFFIX_PARAM].T,linewidth=.6)
            
        if 'MLD_ADJUSTED' in ds.variables:
            ax.scatter(ds[PARAM+SUFFIX_PARAM].sel(N_LEVELS=ds.MLD_ADJUSTED),ds.MLD_ADJUSTED,s=2,color='k')
            
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
        ds = ds.assign_coords(pressure=(name_coord_levels, ds.PRES_INTERP[0,:].data))
    else: # processing = {'','_ADJUSTED'}
        name_coord_levels = "N_LEVELS"
        ds = ds.assign_coords(pressure=(name_coord_levels, ds.PRES[0,:].data))
    
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
def plot_sections(self,PARAM=['TEMP','PSAL','SIG0'],SUFFIX_PARAM='_INTERP',ax=None,namefig=None,figsize=(10, 10),title=None,rolling=0,plot_mld=True,density_threshold=0.03,pmax=1000,**kwargs):


    ds = self
    
    if isinstance(PARAM,str):
        
        if not ax:
            fig, ax = plt.subplots(1,1,figsize=figsize)
        else:
            fig = ax.get_figure()


        if 'cmap' not in kwargs:
            if PARAM=='TEMP':
                kwargs['cmap'] = cmo.thermal
            elif PARAM=='PSAL':
                kwargs['cmap'] = cmo.haline
            else:
                kwargs['cmap'] = cmo.solar
                
        if 'yincrease' not in kwargs:
            kwargs['yincrease'] = False

        if PARAM=='SIG0':
            ds = ds.add_sigma0(SUFFIX_PARAM=SUFFIX_PARAM)
            
        if PARAM+'_INTERP' not in ds.variables:
            ds = ds.add_interp(PARAM)
            
        da = ds[PARAM+'_INTERP']
        if rolling and ds.dims['N_PROF']>rolling :
            da = da.rolling(N_PROF=rolling, center=True, min_periods=1).mean()

        if ds.dims['N_PROF']>10:
            da.T.plot(ax=ax,**kwargs)    

        if plot_mld:
            mld = compute_mld(ds,SUFFIX_PARAM='_ADJUSTED',density_threshold=density_threshold).interpolate_na(dim='N_PROF')
            if rolling and ds.dims['N_PROF']>rolling:
                mld = mld.rolling(N_PROF=rolling,min_periods=min([3,rolling]),center=True).mean()        
            mld.plot(ax=ax)

        if not title:
            number_profiles = np.sum(N_PARAM(ds,PARAM,SUFFIX_PARAM).data>0)
            title = f"{ds.smru_platform_code}, {PARAM+SUFFIX_PARAM}: {number_profiles} profiles"
            
        ax.set_ylim(pmax)
        ax.set_title(title)
        
        # save figure
        if namefig:
            plt.tight_layout()
            plt.savefig(namefig,dpi=300)
            plt.close()
            return
        else:
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

            # save figure
            if namefig:
                plt.tight_layout()
                plt.savefig(namefig,dpi=300,bbox_inches='tight')
                plt.close()
                return
            else:
                plt.tight_layout()
                return fig, ax


        
        
@add_method(xr.Dataset)
def plot_data_tags(self,SUFFIX_PARAM='_ADJUSTED',namefig=None,pmax=1000):

    ds = self
    
    N_LEVELS = 'N_LEVELS'
    if SUFFIX_PARAM=='_INTERP':
        N_LEVELS = 'N_INTERP'
    
    if 'TEMP'+SUFFIX_PARAM not in ds.variables:
        return None, None
    
    if ds['TEMP'+SUFFIX_PARAM].sum(dim=N_LEVELS).sum().data==0:
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

    ds.plot_profiles('TEMP',SUFFIX_PARAM,ax=ax['T'],pmax=pmax)
    ds.plot_profiles('PSAL',SUFFIX_PARAM,ax=ax['S'],pmax=pmax)
    ds.plot_profiles('SIG0',SUFFIX_PARAM,ax=ax['D'],pmax=pmax)
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
        plt.close()
        return
    
    else:
        return fig, ax


##    function to attach a geographic location to a tag  ##
# determine region for each tag
from scipy.interpolate import RegularGridInterpolator
import regionmask

basins = regionmask.defined_regions.ar6.all
label = basins.names
f_region_mask = RegularGridInterpolator((np.arange(-179.5, 180), np.arange(-89.5, 90)), \
                                        basins.mask(np.arange(-179.5, 180), np.arange(-89.5, 90)).transpose().values, \
                                        method='nearest')

def label_regions(ltags):

    # set a new columns called MASK with a regional label
    ltags["MASK"] = f_region_mask(ltags[['LONGITUDE','LATITUDE']].values)
    ltags["MASK"] = ltags.MASK.map(dict(enumerate(label)))
    
    map_regions = {
        'Southern-Ocean':'Southern Ocean',
        'E.Antarctica':'Southern Ocean',
        'W.Antarctica':'Southern Ocean',
        'N.Pacific-Ocean':'North Pacific',
        'C.North-America':'North Pacific', 
        'W.North-America':'North Pacific',
        'N.W.North-America':'North Pacific',
        'N.Central-America':'North Pacific',
        'S.Central-America':'North Pacific',
        'Russian-Arctic':'North Pacific',
        'Arctic-Ocean':'North Atlantic',
        'N.E.North-America':'North Atlantic',
        'E.North-America':'North Atlantic',
        'Greenland/Iceland':'North Atlantic',
        'N.Atlantic-Ocean':'North Atlantic',
        'N.Europe':'North Atlantic',
        'S.E.South-America':'South Atlantic',
        'S.South-America':'South Atlantic',
        'S.Atlantic-Ocean':'South Atlantic',
        'E.Australia':'South Pacific',
        'S.Australia':'South Pacific',
        'New-Zealand':'South Pacific',
        'S.Pacific-Ocean':'South Pacific',
        'Caribbean':'Tropical Atlantic',
        'N.South-America':'Tropical Atlantic',
        'Equatorial.Atlantic-Ocean':'Tropical Atlantic',
        'N.E.South-America':'Tropical Atlantic',
     }
    ltags['MASK'] = ltags.MASK.map(map_regions)
    ltags.loc[(ltags.MASK=='South Pacific')&(ltags.LONGITUDE<0)&(ltags.LONGITUDE>-100),'MASK'] = 'South Atlantic'
    
    return ltags


# build list metadata in MEOP and return as dataframe
@add_method(xr.Dataset)
def list_metadata(self):

    ds=self
    data = {
        'DEPLOYMENT_CODE': ds.deployment_code,
        'SMRU_PLATFORM_CODE': ds.smru_platform_code,
        'CYCLE_NUMBER': ds['CYCLE_NUMBER'].astype(int),
        'JULD': ds['JULD'],
        'LATITUDE': ds['LATITUDE'].round(4),
        'LONGITUDE': ds['LONGITUDE'].round(4),
        'N_TEMP' : N_PARAM(ds,'TEMP'),
        'N_PSAL' : N_PARAM(ds,'PSAL'),
        'N_CHLA' : N_PARAM(ds,'CHLA')}
    df = pd.DataFrame(data)
    df['JULD'] = pd.to_datetime(df.JULD.astype(str))
    ltags = df.groupby('SMRU_PLATFORM_CODE').median()
    mask = label_regions(ltags).MASK
    df = df.merge(mask,on='SMRU_PLATFORM_CODE')

    return df
            


