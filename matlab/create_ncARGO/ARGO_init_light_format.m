% Dimensions:

DATE_TIME = netcdf.defDim(ncid,'DATE_TIME',14);
STRING256 = netcdf.defDim(ncid,'STRING256',256);
STRING64 = netcdf.defDim(ncid,'STRING64',64);
STRING32 = netcdf.defDim(ncid,'STRING32',32);
STRING16 = netcdf.defDim(ncid,'STRING16',16);
STRING8 = netcdf.defDim(ncid,'STRING8',8);
STRING4 = netcdf.defDim(ncid,'STRING4',4);
STRING2 = netcdf.defDim(ncid,'STRING2',2);
N_PROF = netcdf.defDim(ncid,'N_PROF',Nprof);
N_PARAM = netcdf.defDim(ncid,'N_PARAM',3);
N_LEVELS = netcdf.defDim(ncid,'N_LEVELS',Nlevels);
N_CALIB = netcdf.defDim(ncid,'N_CALIB',1);
%N_HISTORY = netcdf.defDim(ncid,'N_HISTORY',netcdf.getConstant('NC_UNLIMITED'));


% Variables and attributes:

varid = netcdf.defVar(ncid,'DATA_TYPE','char',STRING16);
netcdf.putAtt(ncid,varid,'comment','Data type');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'FORMAT_VERSION','char',STRING4);
netcdf.putAtt(ncid,varid,'comment','File format version');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'HANDBOOK_VERSION','char',STRING4);
netcdf.putAtt(ncid,varid,'comment','Data handbook version');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'REFERENCE_DATE_TIME','char',DATE_TIME);
netcdf.putAtt(ncid,varid,'comment','Date of reference for Julian days');
netcdf.putAtt(ncid,varid,'conventions','YYYYMMDDHHMISS');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'DATE_CREATION','char',DATE_TIME);
netcdf.putAtt(ncid,varid,'comment','Date of file creation');
netcdf.putAtt(ncid,varid,'conventions','YYYYMMDDHHMISS');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'DATE_UPDATE','char',DATE_TIME);
netcdf.putAtt(ncid,varid,'long_name','Date of update of this file');
netcdf.putAtt(ncid,varid,'conventions','YYYYMMDDHHMISS');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'PLATFORM_NUMBER','char',[STRING8 N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','Float unique identifier');
netcdf.putAtt(ncid,varid,'conventions','WMO float identifier : A9IIIII');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'PROJECT_NAME','char',[STRING64 N_PROF]);
netcdf.putAtt(ncid,varid,'comment','Name of the project');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'PI_NAME','char',[STRING64 N_PROF]);
netcdf.putAtt(ncid,varid,'comment','Name of the principal investigator');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'CYCLE_NUMBER','long',[N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','Float cycle number');
netcdf.putAtt(ncid,varid,'units','1');
netcdf.putAtt(ncid,varid,'conventions','0..N, 0 : launch cycle (if exists), 1 : first complete cycle');
netcdf.putAtt(ncid,varid,'_FillValue',int32(99999));

varid = netcdf.defVar(ncid,'DIRECTION','char',[N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','Direction of the station profiles');
netcdf.putAtt(ncid,varid,'conventions','A: ascending profiles, D: descending profiles');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'JULD','double',N_PROF);
netcdf.putAtt(ncid,varid,'long_name','Julian day (UTC) of the station relative to REFERENCE_DATE_TIME');
netcdf.putAtt(ncid,varid,'units','days since 1950-01-01 00:00:00 UTC');
netcdf.putAtt(ncid,varid,'calendar','julian');
netcdf.putAtt(ncid,varid,'conventions','Relative julian days with decimal part (as parts of day)');
netcdf.putAtt(ncid,varid,'_FillValue',double(99999));

varid = netcdf.defVar(ncid,'JULD_LOCATION','double',N_PROF);
netcdf.putAtt(ncid,varid,'long_name','Julian day (UTC) of the location relative to REFERENCE_DATE_TIME');
netcdf.putAtt(ncid,varid,'units','days since 1950-01-01 00:00:00 UTC');
netcdf.putAtt(ncid,varid,'calendar','julian');
netcdf.putAtt(ncid,varid,'conventions','Relative julian days with decimal part (as parts of day)');
netcdf.putAtt(ncid,varid,'_FillValue',double(99999));

varid = netcdf.defVar(ncid,'LATITUDE','double',N_PROF);
netcdf.putAtt(ncid,varid,'long_name','Latitude of the station, best estimate');
netcdf.putAtt(ncid,varid,'units','degree_north');
netcdf.putAtt(ncid,varid,'_FillValue',double(99999));
netcdf.putAtt(ncid,varid,'valid_min',-90);
netcdf.putAtt(ncid,varid,'valid_max',90);

varid = netcdf.defVar(ncid,'LONGITUDE','double',N_PROF);
netcdf.putAtt(ncid,varid,'long_name','Longitude of the station, best estimate');
netcdf.putAtt(ncid,varid,'units','degree_east');
netcdf.putAtt(ncid,varid,'_FillValue',double(99999));
netcdf.putAtt(ncid,varid,'valid_min',-180);
netcdf.putAtt(ncid,varid,'valid_max',180);

varid = netcdf.defVar(ncid,'PRES','float',[ N_LEVELS N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','SEA PRESSURE');
netcdf.putAtt(ncid,varid,'units','decibar');
netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
netcdf.putAtt(ncid,varid,'valid_min',0);
netcdf.putAtt(ncid,varid,'valid_max',12000);
netcdf.putAtt(ncid,varid,'comment','In situ measurement, sea surface = 0');

varid = netcdf.defVar(ncid,'PRES_QC','char',[ N_LEVELS N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','quality flag');
netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'TEMP','float',[ N_LEVELS N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','SEA TEMPERATURE IN SITU ITS-90 SCALE');
netcdf.putAtt(ncid,varid,'units','degree_Celsius');
netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
netcdf.putAtt(ncid,varid,'valid_min',-2);
netcdf.putAtt(ncid,varid,'valid_max',40);
netcdf.putAtt(ncid,varid,'comment','In situ measurement');

varid = netcdf.defVar(ncid,'TEMP_QC','char',[ N_LEVELS N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','quality flag');
netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'PSAL','float',[N_LEVELS N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','PRACTICAL SALINITY');
netcdf.putAtt(ncid,varid,'units','1e-3');
netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
netcdf.putAtt(ncid,varid,'valid_min',0);
netcdf.putAtt(ncid,varid,'valid_max',42);
netcdf.putAtt(ncid,varid,'comment','In situ measurement');

varid = netcdf.defVar(ncid,'PSAL_QC','char',[N_LEVELS N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','quality flag');
netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2');
netcdf.putAtt(ncid,varid,'_FillValue',' ');


