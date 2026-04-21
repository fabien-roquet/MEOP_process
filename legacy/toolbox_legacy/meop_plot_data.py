from pathlib import Path
import os
import xarray as xr
import pandas as pd
import numpy as np
import csv
import gsw
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from importlib import reload
import meop
import meop_filenames

plt.rcParams.update({'font.size': 12})

# list functions

# utils
#  central_longitude(ds)

# plot diags
#  plot_data_tags(ds,namefig=None)
#  plot_map_deployments(list_data,namefig=None)
#  plot_data_deployments(list_ds,namefig=None)
#  plot_map_stereo_south(list_ds,namefig=None)

# functions to show matlab-generated plots
#  show_png(figname='')
#  show_plots_depl(EXP,suffix='_hr2')
#  show_plots_tag(one_smru_name,suffix='_hr2')
#  scroll_descriptive_plots_depl(list_files)
#  scroll_calibration_plots_depl(list_files)


#-----------------------------------     utils   --------------------------------------------#
# compute the value of central_longitude in order to create a map, based on the distribution of longitudes
def central_longitude(ds):
    # ds can be a xarray dataset, a list of ds or a dictionary of ds
    if isinstance(ds, list):
        ds2 = []
        for dsi in ds:
            if isinstance(dsi,xr.Dataset):
                ds2.append(dsi['LONGITUDE'])
        lon = xr.concat(ds2, dim="N_PROF")
    elif isinstance(ds, dict):
        ds2 = []
        for tag in ds.keys():
            if isinstance(ds[tag],xr.Dataset):
                ds2.append(ds[tag]['LONGITUDE'])
        lon = xr.concat(ds2, dim="N_PROF")
    else:
        lon = ds['LONGITUDE']
        
    if (lon.max()-lon.min()>180) & (np.abs(lon)>90).all():
        return 180
    else:
        return 0
        
        

def plot_map_deployments(df,groupby='SMRU_PLATFORM_CODE',rebuild=False,
                         namefig=None,folder='.',
                         show_plot=True,
                         title='',legend=True,legend_horiz=False,
                         figsize=(15, 15)
                        ):

    if namefig and (Path(folder) / namefig).exists() and (not rebuild) and (not show_plot):
        return None

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=central_longitude(df)))
    
    df = df.reset_index()
    list_group = df[groupby].unique()
    colors = [
        '#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78', '#2ca02c', '#98df8a',
        '#d62728', '#ff9896', '#9467bd', '#c5b0d5', '#8c564b', '#c49c94',
        '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7', '#bcbd22', '#dbdb8d',
        '#17becf', '#9edae5']
    colors = plt.cm.jet(np.linspace(0, 1, len(list_group)))
    #prop_cycle = plt.rcParams['axes.prop_cycle']
    #colors = prop_cycle.by_key()['color']
    
    dict_cmap = {list_group[i]: colors[i%len(colors)] for i in range(len(list_group))}
    if 'SMRU_PLATFORM_CODE' not in groupby:
        groupby = [groupby,'SMRU_PLATFORM_CODE']
    
    grouped = df.groupby(groupby)
    dict_label = {}
    for name,group in grouped:
        if len(name) == 2:
            name_group = name[0]
        else:
            name_group = name
        if name_group in dict_label:
            label = ''
        else:
            label = name_group
            dict_label[name_group] = 1
        lon = np.where(np.abs(group['LONGITUDE'].diff())>100,np.nan,group['LONGITUDE'])
        lat = group['LATITUDE']
        ax.plot(lon,lat,transform=ccrs.PlateCarree(),color=dict_cmap[name_group],label=label, linewidth=1)

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
              linewidth=2, color='gray', alpha=0.5, linestyle='--')
    #ax.set_extent([-90, -25, -60, -30], ccrs.PlateCarree())
    ax.stock_img()
    ax.coastlines()
    ax.gridlines(xlocs=ax.get_xticks(),ylocs=ax.get_yticks())

    if legend:
        if legend_horiz:
            ncol=len(dict_cmap)
            if ncol < 10:
                ax.legend(bbox_to_anchor=(0, -0.05, 1, 0), loc="upper left", ncol = 10, fontsize=12)
            elif ncol < 50:
                ax.legend(bbox_to_anchor=(0, -0.05, 1, 0), loc="upper left", mode="expand", ncol = 10, fontsize=12)
            else:
                ax.legend(bbox_to_anchor=(0, -0.1, 1, 0), loc="upper left", mode="expand", ncol = 10, fontsize=12)
        else:
            ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left', fontsize=12)

    if title:
        ax.set_title(title, fontsize=14)

    plt.tight_layout()
    if namefig:
        plt.savefig((Path(folder) / namefig),dpi=300,bbox_inches='tight')

    if show_plot:
        return fig, ax
    else:
        plt.close(fig=fig)
        return None



def plot_data_deployments(deployment,namefig=None,rebuild=False,list_fname_prof=None):
    
    if namefig and Path(namefig).exists() and (not rebuild):
        return
    
    if not list_fname_prof:
        list_fname_prof = meop_filenames.list_fname_prof(deployment=deployment,qf='hr1')
    
    from cycler import cycler
    cmap = plt.get_cmap('viridis',len(list_fname_prof))
    custom_cycler = cycler(color=cmap.colors)

    fig = plt.figure(figsize=(12,9))
    gs = fig.add_gridspec(2,3)
    
    ax={}
    ax['T'] = fig.add_subplot(gs[1, 0])
    ax['S'] = fig.add_subplot(gs[1, 1])
    ax['D'] = fig.add_subplot(gs[1, 2])
    ax['TS'] = fig.add_subplot(gs[0, 2:])
    with meop.open_dataset(list_fname_prof[0]) as ds:
        ax['xy'] = fig.add_subplot(gs[0, :2],projection=ccrs.PlateCarree(ds.central_longitude()))
    for key in ax:
        ax[key].set_prop_cycle(custom_cycler)

    SUFFIX_PARAM='_ADJUSTED'
    start_date = []
    end_date = []
    nT, nS, nTS = 0, 0, 0

    list_ds = []
    for fname_prof in list_fname_prof:
        ds = meop.open_dataset(fname_prof)
        ds = ds.add_N_PARAM()
        list_ds.append(ds)
    
    for kk,ds in enumerate(list_ds):
        ds.plot_profiles('TEMP',SUFFIX_PARAM,ax=ax['T'],color=cmap.colors[kk])
        ds.plot_profiles('PSAL',SUFFIX_PARAM,ax=ax['S'],color=cmap.colors[kk])
        ds.plot_profiles('SIG0',SUFFIX_PARAM,ax=ax['D'],color=cmap.colors[kk])
        ds.plot_TSdiag(SUFFIX_PARAM,ax=ax['TS'],draw_sigma=False,color=cmap.colors[kk])
        ds.plot_map(ax=ax['xy'],draw_background=False,color=cmap.colors[kk],scatter_plot=False)
        start_date.append(ds.JULD.min())
        end_date.append(ds.JULD.max())
        nT += np.sum((ds.N_TEMP)>0).data
        nS += np.sum((ds.N_PSAL)>0).data
        nTS+= np.sum((ds.N_TEMP*ds.N_PSAL)>0).data
        ds.close()

    for ds in list_ds:
        ds.close()
        
    ax['T'].set_title(f"TEMP: {nT} T-profiles")
    ax['S'].set_title(f"PSAL: {nS} S-profiles")
    ax['D'].set_title(f"SIG0: {nTS} TS-profiles")

    x0,x1 = ax['TS'].get_xlim()
    y0,y1 = ax['TS'].get_ylim()
    ax['TS'].set_aspect((x1-x0)/(y1-y0))
    X,Y = np.meshgrid(np.linspace(x0,x1),np.linspace(y0,y1))
    D = gsw.sigma0(X,Y)
    ax['TS'].contour(X,Y,D,colors='k')

    gl = ax['xy'].gridlines(crs=ccrs.PlateCarree(), draw_labels=True,
          linewidth=1, color='gray', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    ax['xy'].add_feature(cfeature.LAND, color="lightgrey", alpha=0.5)
    ax['xy'].coastlines()

    # main title
    start_date = min(start_date)
    end_date = min(end_date)
    num_days = (end_date-start_date).dt.days
    ax['xy'].set_title(f"{deployment}: from {start_date.dt.strftime('%Y-%m-%d').data} to {end_date.dt.strftime('%Y-%m-%d').data} ({num_days.data} days)")

    # finitions
    plt.tight_layout()

    # save figure
    if namefig:
        plt.savefig(namefig,dpi=300)
        plt.close(fig=fig)
        return
    
    else:
        return fig, ax



def plot_map_stereo_south(ds,namefig=None,groupby='SMRU_PLATFORM_CODE',title='',legend=True):

    ds = ds.reset_index()
    list_group = sorted(ds[groupby].unique())
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    dict_cmap = {list_group[i]: colors[i%len(colors)] for i in range(len(list_group))}
    if 'SMRU_PLATFORM_CODE' not in groupby:
        groupby = [groupby,'SMRU_PLATFORM_CODE']
    
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.SouthPolarStereo())

    # Compute a circle in axes coordinates, which we can use as a boundary
    # for the map. We can pan/zoom as much as we like - the boundary will be
    # permanently circular.
    import matplotlib.path as mpath
    theta = np.linspace(0, 2*np.pi, 100)
    center, radius = [0.5, 0.5], 0.5
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    circle = mpath.Path(verts * radius + center)
    ax.set_boundary(circle, transform=ax.transAxes)

    grouped = ds.groupby(groupby)
    dict_label = {}
    for name,group in grouped:
        if len(name) == 2:
            name_group = name[0]
        else:
            name_group = name
        if name_group in dict_label:
            label = ''
        else:
            label = name_group
            dict_label[name_group] = 1
        lon = np.where(np.abs(group['LONGITUDE'].diff())>100,np.nan,group['LONGITUDE'])
        lat = group['LATITUDE']
        ax.plot(lon,lat,transform=ccrs.PlateCarree(),color=dict_cmap[name_group],label=label)

    ax.set_extent([-180, 180, -90, -30], ccrs.PlateCarree())
    ax.stock_img()
    ax.coastlines()
    ax.gridlines(xlocs=ax.get_xticks(),ylocs=ax.get_yticks())
    if legend:
        ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    if title:
        ax.set_title(title)

    plt.tight_layout()
    if namefig:
        plt.savefig(namefig,dpi=300,bbox_inches='tight')

    return fig, ax

    

#--------------------------------------------------------------------------------------------#
# functions to show matlab-generated plots
from IPython.display import display,Image
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets
from ipywidgets import HBox

# display a single png figure
def show_png(figname=''):
    figname = Path(figname)
    if figname.exists():
        image = widgets.widgets.Image(value=open(figname, 'rb').read())
        display(HBox([image]))
    else:
        print('File not found:',str(figname))
     
    
# display overview plots for a deployment
def show_plots_depl(deployment,qf='_hr2'):
    mode = ['_recapARGO','_histoARGO']
    figname_mode0 = Path(meop.processdir,'plots',deployment,deployment+qf+'_recapARGO_adj.png')
    figname_mode1 = Path(str(figname_mode0).replace(mode[0], mode[1]))
    if figname_mode0.exists() and figname_mode1.exists():
        imageA = widgets.widgets.Image(value=open(figname_mode0, 'rb').read())
        imageB = widgets.widgets.Image(value=open(figname_mode1, 'rb').read())
        display(HBox([imageA, imageB]))
    return


# display overview plots for a tag
def show_plots_tag(smru_name,qf='_hr2'):
    mode = ['_profiles','_sections']
    EXP = meop.deployment_from_smru_name(smru_name)
    figname_mode0 = Path(meop.processdir,'plots',EXP,smru_name+qf+mode[0]+'.png')
    figname_mode1 = Path(meop.processdir,'plots',EXP,smru_name+qf+mode[1]+'.png')
    if figname_mode0.exists() and figname_mode1.exists():
        imageA = widgets.widgets.Image(value=open(figname_mode0, 'rb').read())
        imageB = widgets.widgets.Image(value=open(figname_mode1, 'rb').read())
        display(HBox([imageA, imageB]))
    return


# show presentation plots for a list of tag from a given deployment
def scroll_descriptive_plots_depl(list_files):

    a = widgets.BoundedIntText(
        value=0,
        min=0,
        max=len(list_files)-1,
        step=1,
        description='Figure:',
        disabled=False
    )

    b = widgets.Dropdown(
        options=[(s.stem,i) for i,s in enumerate(list_files)],
        value=0,
        description='file:',
        disabled=False
    )

    mylink = widgets.link((a, 'value'), (b, 'value'))
    mode = ['_diags_TS','_transect']

    @interact(k=a,k1=b)
    def g(k,k1):
        figname_mode0 = list_files[k]
        figname_mode1 = Path(str(figname_mode0).replace(mode[0], mode[1]))
        if figname_mode0.exists() and figname_mode1.exists():
            imageA = widgets.widgets.Image(value=open(figname_mode0, 'rb').read())
            imageB = widgets.widgets.Image(value=open(figname_mode1, 'rb').read())
            display(HBox([imageA, imageB]))
        return

    
# show calibration plots for a list of tag from a given deployment
def scroll_calibration_plots_depl(list_files):

    a = widgets.BoundedIntText(
        value=0,
        min=0,
        max=len(list_files)-1,
        step=1,
        description='Figure:',
        disabled=False
    )

    b = widgets.Dropdown(
        options=[(str(s.stem).split('_')[1],i) for i,s in enumerate(list_files)],
        value=0,
        description='file:',
        disabled=False
    )

    mylink = widgets.link((a, 'value'), (b, 'value'))
    mode = ['_0','_other_tags']

    @interact(k=a,k1=b)
    def g(k,k1):
        figname_mode0 = list_files[k]
        figname_mode1 = Path(str(figname_mode0).replace(mode[0], mode[1]))
        imageA,imageB=[],[]
        if figname_mode0.exists():
            imageA = widgets.widgets.Image(value=open(figname_mode0, 'rb').read())
            if figname_mode1.exists():
                imageB = widgets.widgets.Image(value=open(figname_mode1, 'rb').read())
                display(HBox([imageA, imageB]))
            else:
                display(HBox([imageA, imageA]))
        return
    
   

