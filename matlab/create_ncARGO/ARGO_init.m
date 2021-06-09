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
N_PARAM = netcdf.defDim(ncid,'N_PARAM',3+isfluo+isoxy+islight);
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

varid = netcdf.defVar(ncid,'STATION_PARAMETERS','char',[STRING16 N_PARAM N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','List of available parameters for the station');
netcdf.putAtt(ncid,varid,'conventions','Argo reference table 3');
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

varid = netcdf.defVar(ncid,'DATA_CENTRE','char',[ STRING2 N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','Data centre in charge of float data processing');
netcdf.putAtt(ncid,varid,'conventions','Argo reference table 4');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'DC_REFERENCE','char',[ STRING32 N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','Station unique identifier in data centre');
netcdf.putAtt(ncid,varid,'conventions','Data centre convention');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'DATA_STATE_INDICATOR','char',[ STRING4 N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','Degree of processing the data have passed through');
netcdf.putAtt(ncid,varid,'conventions','Argo reference table 6');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'DATA_MODE','char',[N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','Delayed mode or real time data');
netcdf.putAtt(ncid,varid,'conventions','R : real time; D : delayed mode; A : real time with adjustment');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'INST_REFERENCE','char',[ STRING64 N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','Instrument type');
netcdf.putAtt(ncid,varid,'conventions','Brand, type, serial number');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'WMO_INST_TYPE','char',[ STRING4 N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','Coded instrument type');
netcdf.putAtt(ncid,varid,'conventions','Argo reference table 8');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'JULD','double',N_PROF);
netcdf.putAtt(ncid,varid,'long_name','Julian day (UTC) of the station relative to REFERENCE_DATE_TIME');
netcdf.putAtt(ncid,varid,'units','days since 1950-01-01 00:00:00 UTC');
netcdf.putAtt(ncid,varid,'calendar','julian');
netcdf.putAtt(ncid,varid,'conventions','Relative julian days with decimal part (as parts of day)');
netcdf.putAtt(ncid,varid,'_FillValue',double(99999));

varid = netcdf.defVar(ncid,'JULD_QC','char',N_PROF);
netcdf.putAtt(ncid,varid,'long_name','Quality on Date and Time');
netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

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

varid = netcdf.defVar(ncid,'POSITION_QC','char',N_PROF);
netcdf.putAtt(ncid,varid,'long_name','Quality on position (latitude and longitude)');
netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'POSITIONING_SYSTEM','char',[ STRING8 N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','Positioning system');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'PROFILE_PRES_QC','char',N_PROF);
netcdf.putAtt(ncid,varid,'long_name','Global quality flag of PRES profile');
netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2a');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'PROFILE_PSAL_QC','char',N_PROF);
netcdf.putAtt(ncid,varid,'long_name','Global quality flag of PSAL profile');
netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2a');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'PROFILE_TEMP_QC','char',N_PROF);
netcdf.putAtt(ncid,varid,'long_name','Global quality flag of TEMP profile');
netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2a');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

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

varid = netcdf.defVar(ncid,'PRES_ADJUSTED','float',[ N_LEVELS N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','SEA PRESSURE');
netcdf.putAtt(ncid,varid,'units','decibar');
netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
netcdf.putAtt(ncid,varid,'valid_min',0);
netcdf.putAtt(ncid,varid,'valid_max',12000);
netcdf.putAtt(ncid,varid,'comment','In situ measurement, sea surface = 0');

varid = netcdf.defVar(ncid,'PRES_ADJUSTED_QC','char',[ N_LEVELS N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','quality flag');
netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'PRES_ADJUSTED_ERROR','float',[ N_LEVELS N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','SEA PRESSURE');
netcdf.putAtt(ncid,varid,'units','decibar');
netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
netcdf.putAtt(ncid,varid,'comment','Contains the error on the adjusted values as determined by the delayed mode QC process.');

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

varid = netcdf.defVar(ncid,'TEMP_ADJUSTED','float',[ N_LEVELS N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','SEA TEMPERATURE IN SITU ITS-90 SCALE');
netcdf.putAtt(ncid,varid,'units','degree_Celsius');
netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
netcdf.putAtt(ncid,varid,'valid_min',-2);
netcdf.putAtt(ncid,varid,'valid_max',40);
netcdf.putAtt(ncid,varid,'comment','In situ measurement');

varid = netcdf.defVar(ncid,'TEMP_ADJUSTED_QC','char',[ N_LEVELS N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','quality flag');
netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'TEMP_ADJUSTED_ERROR','float',[N_LEVELS N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','SEA TEMPERATURE ERROR IN SITU ITS-90 SCALE');
netcdf.putAtt(ncid,varid,'units','degree_Celsius');
netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
netcdf.putAtt(ncid,varid,'comment','Contains the error on the adjusted values as determined by the delayed mode QC process.');

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

varid = netcdf.defVar(ncid,'PSAL_ADJUSTED','float',[ N_LEVELS N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','ADJUSTED PRACTICAL SALINITY');
netcdf.putAtt(ncid,varid,'units','1e-3');
netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
netcdf.putAtt(ncid,varid,'valid_min',0);
netcdf.putAtt(ncid,varid,'valid_max',42);
netcdf.putAtt(ncid,varid,'comment','In situ measurement');

varid = netcdf.defVar(ncid,'PSAL_ADJUSTED_QC','char',[ N_LEVELS N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','quality flag');
netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'PSAL_ADJUSTED_ERROR','float',[ N_LEVELS N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','PRACTICAL SALINITY ERROR');
netcdf.putAtt(ncid,varid,'units','1e-3');
netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
netcdf.putAtt(ncid,varid,'comment','Contains the error on the adjusted values as determined by the delayed mode QC process.');
if issalcor
    varid = netcdf.defVar(ncid,'PSAL_CORRECTED','float',[ N_LEVELS N_PROF]);
    netcdf.putAtt(ncid,varid,'long_name','CORRECTED PRACTICAL SALINITY WITH TIMELAPS');
    netcdf.putAtt(ncid,varid,'units','1e-3');
    netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
    netcdf.putAtt(ncid,varid,'valid_min',0);
    netcdf.putAtt(ncid,varid,'valid_max',42);
    netcdf.putAtt(ncid,varid,'comment','In situ measurement'); 
end
if isfluo
    
    varid = netcdf.defVar(ncid,'PROFILE_CHLA_QC','char',N_PROF);
    netcdf.putAtt(ncid,varid,'long_name','Global quality flag of CHLA profile');
    netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2a');
    netcdf.putAtt(ncid,varid,'_FillValue',' ');
    
    varid = netcdf.defVar(ncid,'CHLA','float',[ N_LEVELS N_PROF]);
    netcdf.putAtt(ncid,varid,'long_name','CHLOROPHYLL-A');
    netcdf.putAtt(ncid,varid,'units','mg/m3');
    netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
    netcdf.putAtt(ncid,varid,'valid_min',0);
    netcdf.putAtt(ncid,varid,'valid_max',10);
    netcdf.putAtt(ncid,varid,'comment','In situ measurement');
    
    varid = netcdf.defVar(ncid,'CHLA_QC','char',[ N_LEVELS N_PROF]);
    netcdf.putAtt(ncid,varid,'long_name','quality flag');
    netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2');
    netcdf.putAtt(ncid,varid,'_FillValue',' ');
    
    varid = netcdf.defVar(ncid,'CHLA_ADJUSTED','float',[ N_LEVELS N_PROF]);
    netcdf.putAtt(ncid,varid,'long_name','CHLOROPHYLL-A');
    netcdf.putAtt(ncid,varid,'units','mg/m3');
    netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
    netcdf.putAtt(ncid,varid,'valid_min',0);
    netcdf.putAtt(ncid,varid,'valid_max',10);
    netcdf.putAtt(ncid,varid,'comment','In situ measurement');
    
    varid = netcdf.defVar(ncid,'CHLA_ADJUSTED_QC','char',[ N_LEVELS N_PROF]);
    netcdf.putAtt(ncid,varid,'long_name','quality flag');
    netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2');
    netcdf.putAtt(ncid,varid,'_FillValue',' ');
    
    varid = netcdf.defVar(ncid,'CHLA_ADJUSTED_ERROR','float',[N_LEVELS N_PROF]);
    netcdf.putAtt(ncid,varid,'long_name','CHLOROPHYLL-A');
    netcdf.putAtt(ncid,varid,'units','mg/m3');
    netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
    netcdf.putAtt(ncid,varid,'comment','Contains the error on the adjusted values as determined by the delayed mode QC process.');

end

if isoxy

    varid = netcdf.defVar(ncid,'PROFILE_DOXY_QC','char',N_PROF);
    netcdf.putAtt(ncid,varid,'long_name','Global quality flag of DOXY profile');
    netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2a');
    netcdf.putAtt(ncid,varid,'_FillValue',' ');
    
    varid = netcdf.defVar(ncid,'DOXY','float',[ N_LEVELS N_PROF]);
    netcdf.putAtt(ncid,varid,'long_name','DISSOLVED OXYGEN');
    netcdf.putAtt(ncid,varid,'units','micromole/kg');
    netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
    netcdf.putAtt(ncid,varid,'valid_min',0);
    netcdf.putAtt(ncid,varid,'valid_max',600);
    netcdf.putAtt(ncid,varid,'comment','In situ measurement');
    
    varid = netcdf.defVar(ncid,'DOXY_QC','char',[ N_LEVELS N_PROF]);
    netcdf.putAtt(ncid,varid,'long_name','quality flag');
    netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2');
    netcdf.putAtt(ncid,varid,'_FillValue',' ');
    
    varid = netcdf.defVar(ncid,'DOXY_ADJUSTED','float',[ N_LEVELS N_PROF]);
    netcdf.putAtt(ncid,varid,'long_name','DISSOLVED OXYGEN');
    netcdf.putAtt(ncid,varid,'units','micromole/kg');
    netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
    netcdf.putAtt(ncid,varid,'valid_min',0);
    netcdf.putAtt(ncid,varid,'valid_max',600);
    netcdf.putAtt(ncid,varid,'comment','In situ measurement');
    
    varid = netcdf.defVar(ncid,'DOXY_ADJUSTED_QC','char',[ N_LEVELS N_PROF]);
    netcdf.putAtt(ncid,varid,'long_name','quality flag');
    netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2');
    netcdf.putAtt(ncid,varid,'_FillValue',' ');
    
    varid = netcdf.defVar(ncid,'DOXY_ADJUSTED_ERROR','float',[N_LEVELS N_PROF]);
    netcdf.putAtt(ncid,varid,'long_name','DISSOLVED OXYGEN');
    netcdf.putAtt(ncid,varid,'units','micromole/kg');
    netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
    netcdf.putAtt(ncid,varid,'comment','Contains the error on the adjusted values as determined by the delayed mode QC process.');
end

if islight

    varid = netcdf.defVar(ncid,'PROFILE_LIGHT_QC','char',N_PROF);
    netcdf.putAtt(ncid,varid,'long_name','Global quality flag of LIGHT profile');
    netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2a');
    netcdf.putAtt(ncid,varid,'_FillValue',' ');
    
    varid = netcdf.defVar(ncid,'LIGHT','float',[ N_LEVELS N_PROF]);
    netcdf.putAtt(ncid,varid,'long_name','ln(PPFD)');
    netcdf.putAtt(ncid,varid,'units','ln(µmol/m²/s)');
    netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
    netcdf.putAtt(ncid,varid,'valid_min',0);
    netcdf.putAtt(ncid,varid,'valid_max',600);
    netcdf.putAtt(ncid,varid,'comment','In situ measurement');
    
    varid = netcdf.defVar(ncid,'LIGHT_QC','char',[ N_LEVELS N_PROF]);
    netcdf.putAtt(ncid,varid,'long_name','quality flag');
    netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2');
    netcdf.putAtt(ncid,varid,'_FillValue',' ');
    
    varid = netcdf.defVar(ncid,'LIGHT_ADJUSTED','float',[ N_LEVELS N_PROF]);
    netcdf.putAtt(ncid,varid,'long_name','ln(PPFD)');
    netcdf.putAtt(ncid,varid,'units','ln(µmol/m²/s)');
    netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
    netcdf.putAtt(ncid,varid,'valid_min',0);
    netcdf.putAtt(ncid,varid,'valid_max',600);
    netcdf.putAtt(ncid,varid,'comment','In situ measurement');
    
    varid = netcdf.defVar(ncid,'LIGHT_ADJUSTED_QC','char',[ N_LEVELS N_PROF]);
    netcdf.putAtt(ncid,varid,'long_name','quality flag');
    netcdf.putAtt(ncid,varid,'conventions','Argo reference table 2');
    netcdf.putAtt(ncid,varid,'_FillValue',' ');
    
    varid = netcdf.defVar(ncid,'LIGHT_ADJUSTED_ERROR','float',[N_LEVELS N_PROF]);
    netcdf.putAtt(ncid,varid,'long_name','ln(PPFD)');
    netcdf.putAtt(ncid,varid,'units','ln(µmol/m²/s)');
    netcdf.putAtt(ncid,varid,'_FillValue',single(99999));
    netcdf.putAtt(ncid,varid,'comment','Contains the error on the adjusted values as determined by the delayed mode QC process.');
end

varid = netcdf.defVar(ncid,'PARAMETER','char',[ STRING16 N_PARAM N_CALIB N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','List of parameters with calibration information');
netcdf.putAtt(ncid,varid,'conventions','Argo reference table 3');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'SCIENTIFIC_CALIB_EQUATION','char',[ STRING256 N_PARAM N_CALIB N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','Calibration equation for this parameter');
netcdf.putAtt(ncid,varid,'_FillValue',' ');

varid = netcdf.defVar(ncid,'SCIENTIFIC_CALIB_COEFFICIENT','char',[ STRING256 N_PARAM N_CALIB N_PROF]);
netcdf.putAtt(ncid,varid,'long_name','Calibration coefficients for this equation');
netcdf.putAtt(ncid,varid,'_FillValue',' ');


% varid = netcdf.defVar(ncid,'SCIENTIFIC_CALIB_COMMENT','char',[ STRING256 N_PARAM N_CALIB N_PROF]);
% netcdf.putAtt(ncid,varid,'long_name','Comment applying to this parameter calibration');
% netcdf.putAtt(ncid,varid,'_FillValue',' ');
%
% varid = netcdf.defVar(ncid,'SCIENTIFIC_CALIB_DATE','char',[ DATE_TIME N_PARAM N_CALIB N_PROF]);
% netcdf.putAtt(ncid,varid,'long_name','Date of calibration');
% netcdf.putAtt(ncid,varid,'_FillValue',' ');



% nc{'HISTORY_INSTITUTION'} = ncchar('N_HISTORY', 'N_PROF', 'STRING4'); %% 0 elements.
% nc{'HISTORY_INSTITUTION'}.long_name = ncchar('Institution which performed action');
% nc{'HISTORY_INSTITUTION'}.conventions = ncchar('Argo reference table 4');
% nc{'HISTORY_INSTITUTION'}._FillValue = ncchar(' ');
%
% nc{'HISTORY_STEP'} = ncchar('N_HISTORY', 'N_PROF', 'STRING4'); %% 0 elements.
% nc{'HISTORY_STEP'}.long_name = ncchar('Step in data processing');
% nc{'HISTORY_STEP'}.conventions = ncchar('Argo reference table 12');
% nc{'HISTORY_STEP'}._FillValue = ncchar(' ');
%
% nc{'HISTORY_SOFTWARE'} = ncchar('N_HISTORY', 'N_PROF', 'STRING4'); %% 0 elements.
% nc{'HISTORY_SOFTWARE'}.long_name = ncchar('Name of software which performed action');
% nc{'HISTORY_SOFTWARE'}.conventions = ncchar('Institution dependent');
% nc{'HISTORY_SOFTWARE'}._FillValue = ncchar(' ');
%
% nc{'HISTORY_SOFTWARE_RELEASE'} = ncchar('N_HISTORY', 'N_PROF', 'STRING4'); %% 0 elements.
% nc{'HISTORY_SOFTWARE_RELEASE'}.long_name = ncchar('Version/release of software which performed action');
% nc{'HISTORY_SOFTWARE_RELEASE'}.conventions = ncchar('Institution dependent');
% nc{'HISTORY_SOFTWARE_RELEASE'}._FillValue = ncchar(' ');
%
% nc{'HISTORY_REFERENCE'} = ncchar('N_HISTORY', 'N_PROF', 'STRING64'); %% 0 elements.
% nc{'HISTORY_REFERENCE'}.long_name = ncchar('Reference of database');
% nc{'HISTORY_REFERENCE'}.conventions = ncchar('Institution dependent');
% nc{'HISTORY_REFERENCE'}._FillValue = ncchar(' ');
%
% nc{'HISTORY_DATE'} = ncchar('N_HISTORY', 'N_PROF', 'DATE_TIME'); %% 0 elements.
% nc{'HISTORY_DATE'}.long_name = ncchar('Date the history record was created');
% nc{'HISTORY_DATE'}.conventions = ncchar('YYYYMMDDHHMISS');
% nc{'HISTORY_DATE'}._FillValue = ncchar(' ');
%
% nc{'HISTORY_ACTION'} = ncchar('N_HISTORY', 'N_PROF', 'STRING4'); %% 0 elements.
% nc{'HISTORY_ACTION'}.long_name = ncchar('Action performed on data');
% nc{'HISTORY_ACTION'}.conventions = ncchar('Argo reference table 7');
% nc{'HISTORY_ACTION'}._FillValue = ncchar(' ');
%
% nc{'HISTORY_PARAMETER'} = ncchar('N_HISTORY', 'N_PROF', 'STRING16'); %% 0 elements.
% nc{'HISTORY_PARAMETER'}.long_name = ncchar('Station parameter action is performed on');
% nc{'HISTORY_PARAMETER'}.conventions = ncchar('Argo reference table 3');
% nc{'HISTORY_PARAMETER'}._FillValue = ncchar(' ');
%
% nc{'HISTORY_START_PRES'} = ncfloat('N_HISTORY', 'N_PROF'); %% 0 elements.
% nc{'HISTORY_START_PRES'}.long_name = ncchar('Start pressure action applied on');
% nc{'HISTORY_START_PRES'}._FillValue = ncfloat(99999);
% nc{'HISTORY_START_PRES'}.units = ncchar('decibar');
%
% nc{'HISTORY_STOP_PRES'} = ncfloat('N_HISTORY', 'N_PROF'); %% 0 elements.
% nc{'HISTORY_STOP_PRES'}.long_name = ncchar('Stop pressure action applied on');
% nc{'HISTORY_STOP_PRES'}._FillValue = ncfloat(99999);
% nc{'HISTORY_STOP_PRES'}.units = ncchar('decibar');
%
% nc{'HISTORY_PREVIOUS_VALUE'} = ncfloat('N_HISTORY', 'N_PROF'); %% 0 elements.
% nc{'HISTORY_PREVIOUS_VALUE'}.long_name = ncchar('Parameter/Flag previous value before action');
% nc{'HISTORY_PREVIOUS_VALUE'}._FillValue = ncfloat(99999);
%
% nc{'HISTORY_QCTEST'} = ncchar('N_HISTORY', 'N_PROF', 'STRING16'); %% 0 elements.
% nc{'HISTORY_QCTEST'}.long_name = ncchar('Documentation of tests performed, tests failed (in hex form)');
% nc{'HISTORY_QCTEST'}.conventions = ncchar('Write tests performed when ACTION=QCP$; tests failed when ACTION=QCF$');
% nc{'HISTORY_QCTEST'}._FillValue = ncchar(' ');
%

