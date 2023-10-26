import meop_process
import meop_metadata
import meop_filenames
import meop
import pandas as pd
from pathlib import Path
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import gsw
import csv

processdir = meop_filenames.processdir


#########################################
print('Process MEOP-CTD')
file_done = processdir / 'python/list_done_viktor.csv'

# create a file to track which tag status
if file_done.is_file():
    ldone = pd.read_csv(file_done,index_col='DEPLOYMENT_CODE')
else:
    lviktor = pd.read_csv(processdir / 'python/list_deployment_viktor.csv',usecols=[1])
    ldone = lviktor.copy()
    ldone['DONE'] = 0
    ldone = ldone.set_index('DEPLOYMENT_CODE')
    ldone.to_csv(file_done)


    
    
#########################################
# process tags
meop_process.start_matlab()
conf = meop_process.init_mirounga()

#ldeployment = meop_metadata.read_list_deployment()
#ldeployment = ldeployment[ldeployment.PROCESS==1]

for deployment in ldone[ldone.DONE.isin([0,9])].index:

    print("Process ",deployment)
    success = True
    #meop_process.import_raw_data(deployment=deployment)
    success *= meop_process.process_tags(deployment=deployment,smru_name='',notlc=True)
    success *= meop_process.create_hr2(deployment=deployment,smru_name='')
    if success:
        ldone.loc[deployment,'DONE']=1
    else:
        ldone.loc[deployment,'DONE']=9
    ldone.to_csv(file_done)

meop_process.stop_matlab()



#########################################

# extract information
def profile_data(ds):
    
    # ds: xarray dataset containing data from a tag
    ds['profileID']=(('N_PROF'),  ['{}_{:04d}'.format(ds.smru_platform_code,kk) for kk in range(ds.dims['N_PROF'])])
    ds['YEAR'] = (('N_PROF'),[0]*ds.dims['N_PROF'])
    ds['MONTH'] = (('N_PROF'),[0]*ds.dims['N_PROF'])
    ds['DAY'] = (('N_PROF'),[0]*ds.dims['N_PROF'])
    ds['HOUR'] = (('N_PROF'),[0]*ds.dims['N_PROF'])
    ds['MINUTE'] = (('N_PROF'),[0]*ds.dims['N_PROF'])
    ds['SECOND'] = (('N_PROF'),[0]*ds.dims['N_PROF'])
    for kk,time in enumerate(ds['JULD'].data):
        ds['YEAR'][kk]=time.year
        ds['MONTH'][kk]=time.month
        ds['DAY'][kk]=time.day
        ds['HOUR'][kk]=time.hour
        ds['MINUTE'][kk]=time.minute
        ds['SECOND'][kk]=time.second
    
    fieldnames = ['profileID','LATITUDE','LONGITUDE','YEAR','MONTH','DAY','HOUR','PRES_ADJUSTED','TEMP_ADJUSTED','PSAL_ADJUSTED','TEMP_QC','PSAL_QC']
    df = ds[fieldnames].to_dataframe()
    df = df.reset_index()
    df = df.set_index(['N_PROF', 'N_LEVELS']).sort_index()

    # remove rows with no data
    I = np.logical_or(df['TEMP_QC']==b'1', df['PSAL_QC']==b'1')
    df = df[I]

    # remove rows with no PRES data
    df = df[~np.isnan(df['PRES_ADJUSTED'])]
    
    # replace nan with 99.999
    df = df.fillna(99.999)
    
    df['PRES_ADJUSTED']=df['PRES_ADJUSTED'].astype(int)
    df['TEMP_QC']=df['TEMP_QC'].astype(int)
    df['PSAL_QC']=df['PSAL_QC'].astype(int)
    
    return df



#########################################
# create csv files for viktor
for deployment in ldone[ldone.DONE.isin([1])].index:

    print("Process ",deployment)
    folder_out = processdir / 'MEOP_csv_gz'
    folder_out.mkdir(parents=True, exist_ok=True)
    
    qf='hr1'
    fnames = meop_filenames.list_fname_prof(deployment=deployment,qf=qf)
    success = True
    for fname in fnames:
        smru_name = meop_filenames.smru_name_from_fname_prof(fname)
        print('...   ',smru_name)        
        
        try:
            ds = meop.open_dataset(fname)
            df = profile_data(ds)
            cols = ['profileID','LATITUDE','LONGITUDE','YEAR','MONTH','DAY','HOUR']
            df.loc[df.duplicated(subset='profileID',keep='first'),cols]=np.nan
            filename = folder_out / (smru_name+'_'+qf+'.csv.gz')
            if not filename.is_file():
                df.to_csv(filename, mode='w', index=False, float_format='%6.3f',compression='gzip')

        except:
            print('         ERROR')
            success = False
        
    if success:
        ldone.loc[deployment,'DONE']=2
        ldone.to_csv(file_done)

    

#########################################
# generate plots
meop_process.start_matlab()
conf = meop_process.init_mirounga()

for deployment in ldone[ldone.DONE.isin([2])].index:

    print("Process ",deployment)
    success = True
    #meop_process.generate_calibration_plots(deployment=deployment,smru_name='')
    meop_process.generate_doc_latex(deployment=deployment,smru_name='')
    #meop_process.export_odv4(deployment=deployment,smru_name='')
    if success:
        ldone.loc[deployment,'DONE']=3
        ldone.to_csv(file_done)

meop_process.stop_matlab()



