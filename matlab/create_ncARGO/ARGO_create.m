function ARGO_create(ficout,Nprof,Nlevels,isfluo,isoxy,islight,issalcor)
%% create a file in ARGO netcdf format


ncid = netcdf.create(ficout, 'clobber');

ARGO_init;

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
a=ncread(ficout,'STATION_PARAMETERS');
for kk=1:Nprof, 
    a(:,1,kk)=sprintf('%-16s','PRES'); 
    a(:,2,kk)=sprintf('%-16s','TEMP'); 
    a(:,3,kk)=sprintf('%-16s','PSAL'); 
    if isfluo, a(:,4,kk)=sprintf('%-16s','CHLA'); end
    if isoxy , a(:,4+isfluo,kk)=sprintf('%-16s','DOXY'); end
    if islight , a(:,4+isfluo+isoxy,kk)=sprintf('%-16s','LIGHT'); end

end
ncwrite(ficout,'STATION_PARAMETERS',a);
ncwrite(ficout,'CYCLE_NUMBER',zeros(Nprof,1)*NaN);
ncwrite(ficout,'DIRECTION',repmat('A',1,Nprof));
ncwrite(ficout,'DATA_CENTRE',repmat('IF',Nprof,1)');
ncwrite(ficout,'DC_REFERENCE',repmat(' ',32,Nprof));
ncwrite(ficout,'DATA_MODE',repmat('D',1,Nprof));
ncwrite(ficout,'WMO_INST_TYPE',repmat('995 ',Nprof,1)');
ncwrite(ficout,'JULD',zeros(Nprof,1)*NaN);
ncwrite(ficout,'JULD_LOCATION',zeros(Nprof,1)*NaN);
ncwrite(ficout,'JULD_QC',repmat('1',1,Nprof));
ncwrite(ficout,'LATITUDE',zeros(Nprof,1)*NaN);
ncwrite(ficout,'LONGITUDE',zeros(Nprof,1)*NaN);
ncwrite(ficout,'POSITION_QC',repmat('1',1,Nprof));
%%
ncwrite(ficout,'POSITIONING_SYSTEM',repmat('ARGOS   ',Nprof,1)');
ncwrite(ficout,'PROFILE_PRES_QC',repmat('A',1,Nprof));
ncwrite(ficout,'PROFILE_PSAL_QC',repmat('A',1,Nprof));
ncwrite(ficout,'PROFILE_TEMP_QC',repmat('A',1,Nprof));
if isfluo, ncwrite(ficout,'PROFILE_CHLA_QC',repmat('A',1,Nprof)); end
if isoxy,  ncwrite(ficout,'PROFILE_DOXY_QC',repmat('A',1,Nprof)); end
if islight,  ncwrite(ficout,'PROFILE_LIGHT_QC',repmat('A',1,Nprof)); end

%%
PARAMETER = ncread(ficout,'PARAMETER');
SCIENTIFIC_CALIB_EQUATION = ncread(ficout,'SCIENTIFIC_CALIB_EQUATION');
SCIENTIFIC_CALIB_COEFFICIENT = ncread(ficout,'SCIENTIFIC_CALIB_COEFFICIENT');
for kk=1:Nprof,
    PARAMETER(:,1,1,kk)=sprintf('%-16s','PRES');
    PARAMETER(:,2,1,kk)=sprintf('%-16s','TEMP');
    PARAMETER(:,3,1,kk)=sprintf('%-16s','PSAL');
    SCIENTIFIC_CALIB_EQUATION(:,1,1,kk)=   sprintf('%-256s','Pc = P - ( p1 [dbar/km] * P  * 1e-3 + p2 [dbar] )');
    SCIENTIFIC_CALIB_EQUATION(:,2,1,kk)=   sprintf('%-256s','Tc = T - ( t1 [degC/km] * Pc * 1e-3 + t2 [degC] )');
    SCIENTIFIC_CALIB_EQUATION(:,3,1,kk)=   sprintf('%-256s','Sc = S - ( s1 [ psu/km] * Pc * 1e-3 + s2 [ psu] )');
    SCIENTIFIC_CALIB_COEFFICIENT(:,1,1,kk)=sprintf('%-256s',' ');
    SCIENTIFIC_CALIB_COEFFICIENT(:,2,1,kk)=sprintf('%-256s',' ');
    SCIENTIFIC_CALIB_COEFFICIENT(:,3,1,kk)=sprintf('%-256s',' ');
    if isfluo
        PARAMETER(:,4,1,kk)=sprintf('%-16s','CHLA');
        SCIENTIFIC_CALIB_EQUATION(:,4,1,kk)=   sprintf('%-256s','Fc = f1 * F + f2');
        SCIENTIFIC_CALIB_COEFFICIENT(:,4,1,kk)=sprintf('%-256s',' ');
    end
    if isoxy
        PARAMETER(:,4+isfluo,1,kk)=sprintf('%-16s','DOXY');
        SCIENTIFIC_CALIB_EQUATION(:,4+isfluo,1,kk)=   sprintf('%-256s',' ');
        SCIENTIFIC_CALIB_COEFFICIENT(:,4+isfluo,1,kk)=sprintf('%-256s',' ');
    end
    if islight
        PARAMETER(:,4+isfluo+isoxy,1,kk)=sprintf('%-16s','LIGHT');
        SCIENTIFIC_CALIB_EQUATION(:,4+isfluo+isoxy,1,kk)=   sprintf('%-256s',' ');
        SCIENTIFIC_CALIB_COEFFICIENT(:,4+isfluo+isoxy,1,kk)=sprintf('%-256s',' ');
    end
end
ncwrite(ficout,'PARAMETER', PARAMETER);
ncwrite(ficout,'SCIENTIFIC_CALIB_EQUATION', SCIENTIFIC_CALIB_EQUATION);
ncwrite(ficout,'SCIENTIFIC_CALIB_COEFFICIENT', SCIENTIFIC_CALIB_COEFFICIENT);




