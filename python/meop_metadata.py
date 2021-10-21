from pathlib import Path
import os
import xarray as xr
import pandas as pd
import numpy as np
import gsw
import netCDF4 as nc
import meop
import meop_filenames

processdir = meop_filenames.processdir

# list functions

#  read_list_profiles(rebuild=False,verbose=False,public=False,Tdata=False,country=None,qf='lr0')

#  label_regions(ltags):
#  filter_public_data(lprofiles, ltags, ldeployments)
#  filter_profiles_with_Tdata(lprofiles, ltags, ldeployments)
#  filter_country(country, lprofiles, ltags, ldeployments)


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



# list_deployment obtained with read_list_deployment()
def list_profiles_from_ncfile(qf='lr0',datadir=(processdir / 'final_dataset_prof'), \
                        deployment='',country='',deployments=[],countries=[], \
                        public_only=False,file_pkl=''):
    
    df = read_list_deployment().reset_index()
    if deployments:
        df = df[df.DEPLOYMENT_CODE.isin(deployments)]
    if deployment:
        df = df[df.DEPLOYMENT_CODE.isin([deployment])]
    if countries:
        df = df[df.COUNTRY.isin(countries)]
    if country:
        df = df[df.COUNTRY.isin([country])]
    if public_only:
        df = df[df.PUBLIC==1]
    
    # list of deployment dataframes
    list_df=[]

    for deployment in df.DEPLOYMENT_CODE:
        list_fname = meop_filenames.list_fname_prof(deployment=deployment,qf=qf)
        
        # concatenate df by deployments
        for ncfile in list_fname:
            if ncfile.exists():
                with meop.open_dataset(ncfile) as ds:
                    df = ds.list_metadata()
                list_df.append(df)
                
    if not list_df:
        return []
    
    # concatenate list of dataframes into one dataframe
    lprofiles = pd.concat(list_df)

    # add label regions
    df_tag = lprofiles.groupby('SMRU_PLATFORM_CODE').mean()
    mask = label_regions(df_tag).MASK
    lprofiles = lprofiles.merge(mask,on='SMRU_PLATFORM_CODE')

    if file_pkl:
        lprofiles.to_pickle(processdir / file_pkl)

    ltags, ldeployments = list_tags_deployments_from_profiles(lprofiles)
    lprofiles[['year','month','day']] = [[lprofiles.JULD[kk].year,lprofiles.JULD[kk].month,lprofiles.JULD[kk].day] for kk in range(len(lprofiles))]
    
    return lprofiles, ltags, ldeployments
   
    
    
# read table files and append information about tags and deployments
def list_tags_deployments_from_profiles(lprofiles):
    
    ltags = lprofiles.groupby('SMRU_PLATFORM_CODE').first()\
        .drop(['N_TEMP','N_PSAL','N_CHLA','CYCLE_NUMBER'],axis='columns')
    ltags['N_PROF_TEMP'] = lprofiles.groupby('SMRU_PLATFORM_CODE').N_TEMP.count()
    ltags['N_PROF_PSAL'] = lprofiles.groupby('SMRU_PLATFORM_CODE').N_PSAL.count()
    ltags['N_PROF_CHLA'] = lprofiles.groupby('SMRU_PLATFORM_CODE').N_CHLA.count()
    
    agg_ops = {'JULD': min, 'LATITUDE': np.mean, 'LONGITUDE': np.mean, 'N_PROF_TEMP': sum,
       'N_PROF_PSAL': sum, 'N_PROF_CHLA': sum}
    ldeployments = ltags.groupby('DEPLOYMENT_CODE').agg(agg_ops)
    ldeployments['N_TAGS'] = ltags.groupby('DEPLOYMENT_CODE').DEPLOYMENT_CODE.count()
    ldeployments = ldeployments.merge(read_list_deployment(),on='DEPLOYMENT_CODE',how='left')
    drop_list = ['START_DATE','END_DATE','START_DATE_JUL']
    ldeployments = ldeployments.drop(drop_list,axis='columns')
    
    list_public = ldeployments.reset_index()[['DEPLOYMENT_CODE','PUBLIC']]
    ltags = ltags.reset_index().merge(list_public,on='DEPLOYMENT_CODE')

    # add correction coefficients in ltags
    if (processdir / 'table_coeff.csv').is_file():
        coeff = pd.read_csv(processdir / 'table_coeff.csv')
        coeff = coeff.rename(columns={'comment':'comments'})
        ltags = ltags.merge(coeff,left_on='SMRU_PLATFORM_CODE',right_on='smru_platform_code',how='left')
        ltags = ltags.drop('smru_platform_code', axis='columns')   
    
    
    # add variable_offset in ltags
    if (processdir / 'table_salinity_offsets.csv').is_file():
        salinity_offsets = pd.read_csv(processdir / 'table_salinity_offsets.csv')
        salinity_offsets['variable_offset'] = 1
        variable_offset = salinity_offsets[['smru_platform_code','variable_offset']]
        ltags = ltags.merge(variable_offset,left_on='SMRU_PLATFORM_CODE',right_on='smru_platform_code',how='left')\
            .drop('smru_platform_code', axis='columns')

    list_deployment_hr = read_list_deployment_hr()
    if not isinstance(list_deployment_hr,list):
        ltags = ltags.merge(list_deployment_hr,left_on='SMRU_PLATFORM_CODE',right_on='smru_platform_code',how='left')        
        ltags = ltags.drop('smru_platform_code', axis='columns')
    
    ltags['comments'] = ltags['comments'].fillna('no comment')
    
    return ltags, ldeployments
    
    
# determine region for each tag
def label_regions(ltags):

    # set a new columns called MASK with a regional label
    from scipy.interpolate import RegularGridInterpolator
    import regionmask

    basins = regionmask.defined_regions.ar6.all
    label = basins.names
    lon = np.arange(-179.5, 180)
    lat = np.arange(-89.5, 90)
    mask = basins.mask(lon,lat)
    f = RegularGridInterpolator((lon, lat), mask.transpose().values,method='nearest')
    ltags["MASK"] = f(ltags[['LONGITUDE','LATITUDE']].values)
    ltags["MASK"] = ltags.MASK.map(dict(enumerate(label)))
    
    map_regions = {
        'Southern-Ocean':'Southern Ocean',
        'E.Antarctica':'Southern Ocean',
        'W.Antarctica':'Southern Ocean',
        'Arctic-Ocean':'North Atlantic',
        'N.Pacific-Ocean':'North Pacific',
        'C.North-America':'North Pacific', 
        'W.North-America':'North Pacific',
        'N.E.North-America':'North Atlantic',
        'E.North-America':'North Atlantic',
        'Greenland/Iceland':'North Atlantic',
        'N.Atlantic-Ocean':'North Atlantic',
        'N.W.North-America':'North Pacific',
        'N.Europe':'North Atlantic',
        'S.Australia':'Australia',
        'N.Central-America':'North Pacific',
        'N.South-America':'Tropical Atlantic',
        'S.South-America':'South Atlantic',
        'S.Atlantic-Ocean':'South Atlantic',
    }
    ltags['MASK'] = ltags.MASK.map(map_regions)

    return ltags


# select only public data
def filter_public_data(lprofiles, ltags, ldeployments):    
    ltags = ltags[ltags.PUBLIC == 1]
    ldeployments = ldeployments[ldeployments.PUBLIC == 1]
    lprofiles = lprofiles.merge(ltags.SMRU_PLATFORM_CODE,on='SMRU_PLATFORM_CODE')
    return lprofiles, ltags, ldeployments


# select only profiles with data points
def filter_profiles_with_Tdata(lprofiles, ltags, ldeployments):    
    ltags = ltags[ltags.N_PROF_TEMP!=0]
    ldeployments = ldeployments.merge(ltags.DEPLOYMENT_CODE,on='DEPLOYMENT_CODE')
    lprofiles = lprofiles.loc[lprofiles.N_TEMP!=0]
    return lprofiles, ltags, ldeployments


# select only profiles with data points
def filter_country(country, lprofiles, ltags, ldeployments):    
    ldeployments = ldeployments.loc[ldeployments.COUNTRY==country]
    ltags = ltags[ltags.DEPLOYMENT_CODE.isin(ldeployments.DEPLOYMENT_CODE)]
    lprofiles = lprofiles[lprofiles.SMRU_PLATFORM_CODE.isin(ltags.SMRU_PLATFORM_CODE)]
    return lprofiles, ltags, ldeployments


# read MEOP data list from pickle file and return the dataframe.
# If filename_pkl is not found, the list file is generated.
def read_lists_metadata(file_pkl='',rebuild=False,\
                       save_to_pkl=True,\
                       verbose=False,public=False,Tdata=False,country=None,
                       qf='lr0',datadir=(processdir / 'final_dataset_prof')):

    if not file_pkl:
        file_pkl = processdir / f'list_profiles_{qf}.pkl'
        
    if file_pkl.is_file() and (not rebuild):
        lprofiles = pd.read_pickle(file_pkl)
        ltags, ldeployments = list_tags_deployments_from_profiles(lprofiles)
    else:
        lprofiles, ltags, ldeployments = list_profiles_from_ncfile(qf=qf, datadir=datadir, public_only=public)
        if save_to_pkl:
            lprofiles.to_pickle(file_pkl)
    
    # track tags with issues
    tag_problem = ltags.loc[ltags.SMRU_PLATFORM_CODE.isnull(),:]
    if len(tag_problem.instr_id):
        print('List of instr id for tags with hr datasets but no low resolution ones:')
        print(tag_problem.instr_id)
        for tag in list(tag_problem.index):
            ltags.drop(tag,axis=0, inplace=True)
        print('')

    tag_problem = ltags.loc[ltags.SMRU_PLATFORM_CODE.isnull(),:]
    if len(tag_problem):
        print('List of tags with correction coefficients yet not listed in list_deployment:')
        print(tag_problem.smru_platform_code)
        message = 'tag with correction coefficient, yet no netcdf file'
        for tag in list(tag_problem.SMRU_PLATFORM_CODE):
            comment = coeff.loc[coeff.smru_platform_code == tag,'comment']            
            if message not in comment:
                if 'no comment' in comment:
                    comment = message
                else:
                    comment = comment+', '+message
                coeff.loc[coeff.smru_platform_code == tag,'comment'] = comment
                ltags.loc[ltags.SMRU_PLATFORM_CODE == tag,'comment'] = comment    
                
    if public:
        lprofiles, ltags, ldeployments = filter_public_data(lprofiles, ltags, ldeployments)

    if Tdata:
        lprofiles, ltags, ldeployments = filter_profiles_with_Tdata(lprofiles, ltags, ldeployments)
        
    if country:
        lprofiles, ltags, ldeployments = filter_country(country, lprofiles, ltags, ldeployments)
        
    return lprofiles, ltags, ldeployments





# Execute in terminal command line
if __name__ == "__main__":

    import argparse
    
    # create parser
    parser = argparse.ArgumentParser()

    # add arguments to the parser
    parser.add_argument("--rebuild", help = "Rebuild list of metadata from ncfiles? [default False]", action='store_true')
    parser.add_argument("--verbose", help = "Verbose mode [default False]", action='store_true')
    parser.add_argument("--public", help = "Load public data only [default False]", action='store_true')
    parser.add_argument("--Tdata", help = "Load only profiles with at least one temperature profile [default False]", action='store_true')
    parser.add_argument("--qf", default='lr0', help = "quality flag: lr0, lr1, hr0, hr1, fr0, fr1, hr2 [default lr0]")

    # parse the arguments
    args = parser.parse_args()

    lprofiles, ltags, ldeployments = \
        read_list_profiles(rebuild=args.rebuild,verbose=args.verbose,public=args.public,Tdata=args.Tdata,qf=args.qf)
    
    

