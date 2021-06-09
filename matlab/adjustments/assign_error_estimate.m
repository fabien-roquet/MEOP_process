function assign_error_estimate(name_prof,name_lr0,temp_error,psal_error)
% error taken from file table_param.csv
% error is multiplied by for low resolution profiles

% set first guess of observational error
error = ncload_struct(name_prof,'TEMP_ADJUSTED_ERROR','PSAL_ADJUSTED_ERROR');
TEMP_QC = double(ncread(name_prof,'TEMP_QC'));
TEMP_QC(TEMP_QC>47) = TEMP_QC(TEMP_QC>47)-48; 
TEMP_QC(TEMP_QC==32) = 9;

PSAL_QC = double(ncread(name_prof,'PSAL_QC'));
PSAL_QC(PSAL_QC>47) = PSAL_QC(PSAL_QC>47)-48; 
PSAL_QC(PSAL_QC==32) = 9;

error.TEMP_ADJUSTED_ERROR(TEMP_QC<9)=temp_error;
error.PSAL_ADJUSTED_ERROR(PSAL_QC<9)=psal_error;
error.TEMP_ADJUSTED_ERROR(TEMP_QC==9)=NaN;
error.PSAL_ADJUSTED_ERROR(PSAL_QC==9)=NaN;

TEMP_QC_LR0 = double(ncread(name_lr0,'TEMP_QC'));
TEMP_QC_LR0(TEMP_QC_LR0>47) = TEMP_QC_LR0(TEMP_QC_LR0>47)-48; 
TEMP_QC_LR0(TEMP_QC_LR0==32) = 9;

PSAL_QC_LR0 = double(ncread(name_lr0,'PSAL_QC'));
PSAL_QC_LR0(PSAL_QC_LR0>47) = PSAL_QC_LR0(PSAL_QC_LR0>47)-48; 
PSAL_QC_LR0(PSAL_QC_LR0==32) = 9;

if TEMP_QC_LR0 & ...
  size(TEMP_QC_LR0,2)==size(error.TEMP_ADJUSTED_ERROR,2) & ...
  size(PSAL_QC,2)==size(error.PSAL_ADJUSTED_ERROR,2),    
  
    nT=nansum(double(TEMP_QC_LR0<=1));
    J = find(nT<10);
    if length(J), 
        error.TEMP_ADJUSTED_ERROR(:,J)=temp_error*2; 
    end
    error.TEMP_ADJUSTED_ERROR(TEMP_QC==9)=NaN;

    nS=nansum(double(PSAL_QC_LR0<=1));
    J = find(nS<10);
    if length(J), 
        error.PSAL_ADJUSTED_ERROR(:,J)=psal_error*2; 
    end
    error.PSAL_ADJUSTED_ERROR(PSAL_QC==9)=NaN;

end

ncwrite(name_prof,'TEMP_ADJUSTED_ERROR',error.TEMP_ADJUSTED_ERROR);
ncwrite(name_prof,'PSAL_ADJUSTED_ERROR',error.PSAL_ADJUSTED_ERROR);


