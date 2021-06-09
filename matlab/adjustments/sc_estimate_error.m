function assign_error_estimate(name_prof,name_lr0,temp_error,psal_error)
% error taken from file table_param.csv
% error is multiplied by for low resolution profiles

Mqc_lr0 = ARGO_load_qc(name_lr0,0);
I=1:length(Mqc_lr0.JULD);%find(strcmp(Mqc.platform_number,Mqc.list_descr{kk}));
if isempty(I), return; end
nT=nansum(double(Mqc_lr0.TEMP_QC(:,:)<=1));
nS=nansum(double(Mqc_lr0.PSAL_QC(:,:)<=1));

% set first guess of observational error
error = ncload_struct(name_prof,'TEMP_ADJUSTED_ERROR','PSAL_ADJUSTED_ERROR');

error.TEMP_ADJUSTED_ERROR(:)=temp_error;
J = find(nT<10);
if length(J), 
    error.TEMP_ADJUSTED_ERROR(:,J)=temp_error*2; 
end
error.TEMP_ADJUSTED_ERROR(Mqc.TEMP_QC==9)=NaN;

error.PSAL_ADJUSTED_ERROR(:)=psal_error;
J = find(nS<10);
if length(J), 
    error.PSAL_ADJUSTED_ERROR(:,J)=psal_error*2; 
end
error.PSAL_ADJUSTED_ERROR(Mqc.PSAL_QC==9)=NaN;

ncwrite(name_prof,'TEMP_ADJUSTED_ERROR',error.TEMP_ADJUSTED_ERROR);
ncwrite(name_prof,'PSAL_ADJUSTED_ERROR',error.PSAL_ADJUSTED_ERROR);


