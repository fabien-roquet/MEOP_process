function ARGO_create(ficout,Nprof,Nlevels)
%% create a file in ARGO netcdf format


ncid = netcdf.create(ficout, 'clobber');

ARGO_init_light_format;

%varid = netcdf.defVar(ncid,'SMRU_NAME','char',[STRING32 N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','SRDL identifier in the SMRU database');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

netcdf.close(ncid)

%% load data:
ncwrite(ficout,'DATA_TYPE','Argo profile    ');
ncwrite(ficout,'FORMAT_VERSION','3.0 ');
ncwrite(ficout,'HANDBOOK_VERSION','3.0 ');
ncwrite(ficout,'REFERENCE_DATE_TIME','19500101000000');
ncwrite(ficout,'DATE_CREATION',datestr(now,'yyyymmddHHMMSS'));
ncwrite(ficout,'DATE_UPDATE',datestr(now,'yyyymmddHHMMSS'));
ncwrite(ficout,'PLATFORM_NUMBER',repmat(' ',8,Nprof));
ncwrite(ficout,'PROJECT_NAME',repmat(sprintf('%64s',' '),Nprof,1)');
ncwrite(ficout,'PI_NAME',repmat(sprintf('%64s',' '),Nprof,1)');
ncwrite(ficout,'CYCLE_NUMBER',zeros(Nprof,1)*NaN);
ncwrite(ficout,'JULD',zeros(Nprof,1)*NaN);
ncwrite(ficout,'JULD_LOCATION',zeros(Nprof,1)*NaN);
ncwrite(ficout,'LATITUDE',zeros(Nprof,1)*NaN);
ncwrite(ficout,'LONGITUDE',zeros(Nprof,1)*NaN);

%%
for kk=1:Nprof,
    PARAMETER(:,1,1,kk)=sprintf('%s%12s','PRES',' ');
    PARAMETER(:,2,1,kk)=sprintf('%s%12s','TEMP',' ');
    PARAMETER(:,3,1,kk)=sprintf('%s%12s','PSAL',' ');
end




