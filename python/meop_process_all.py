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
file_done = processdir / 'python/list_done_all.csv'

# create a file to track which tag status
if file_done.is_file():
    ldone = pd.read_csv(file_done,index_col='DEPLOYMENT_CODE')
else:    
    list_deployment = meop_metadata.read_list_deployment()
    list_deployment = list_deployment[list_deployment.PROCESS==1]
    ldone = list_deployment.reset_index()[['DEPLOYMENT_CODE']]
    ldone['DONE1'] = 0
    ldone['DONE2'] = 0
    ldone['DONE3'] = 0
    ldone = ldone.set_index('DEPLOYMENT_CODE')
    ldone.to_csv(file_done)

    
#########################################
# process tags
meop_process.start_matlab()
conf = meop_process.init_mirounga()

for deployment in ldone[ldone.DONE1.isin([0,9])].index:

    print("Process ",deployment)
    meop_process.import_raw_data(deployment=deployment)
    if (   
       meop_process.process_tags(deployment=deployment) and
       meop_process.create_hr2(deployment=deployment) and
       meop_process.generate_doc_latex(deployment=deployment)
    ):
        ldone.loc[deployment,'DONE1']=1
    else:
        ldone.loc[deployment,'DONE1']=9
    ldone.to_csv(file_done)


for deployment in ldone[ldone.DONE1==1].index:

    if ldone[ldone.DONE2.isin([0,9])].index:
        print("Process ",deployment)
        if (   
           meop_process.generate_calibration_plots(deployment=deployment) and
           meop_process.export_odv4(deployment=deployment)
        ):
            ldone.loc[deployment,'DONE2']=1
        else:
            ldone.loc[deployment,'DONE2']=9
        ldone.to_csv(file_done)

    
meop_process.stop_matlab()


    



