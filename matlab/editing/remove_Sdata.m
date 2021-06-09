%% remove S data

Mqc.PSAL_QC(I)=4;
Mqc.PSAL(Mqc.PSAL_QC>1)=NaN;

nS=nansum(double(Mqc.PSAL_QC<2));
Mqc=remove_Sprofiles(info_deployment,smru_name,'index',find(nS<=4),suffix);

Mqc.PRES_QC(Mqc.PSAL_QC(:)==4 & Mqc.TEMP_QC(:)==9)=4;
Mqc.PRES(Mqc.PRES_QC>1)=NaN;

