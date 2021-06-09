
%% Creation du fichier nc timeseries
ficout = sprintf('%s%s_traj.nc',info_deployment.dir,smru_name);
if exist(ficout,'file'),
    delete(ficout);
end

    
%     prof_pres_cor=P-(P*coef.p1+coef.p2);
%     prof_temp_cor=T-(prof_pres_cor*coef.t1+coef.t2);
%     prof_sal_cor=S-(prof_pres_cor*coef.s1+coef.s2);
Nhr = length(hrdata.P);
prof_pres_cor=hrdata.P;
prof_temp_cor=hrdata.T;
prof_sal_cor=hrdata.S;
prof_fluo_cor=hrdata.F;
prof_oxy_cor=hrdata.O;
latinterp=interp1(jul_hr,lat_hr,hrdata.date,'linear','extrap');
loninterp=interp1(jul_hr,lon_hr,hrdata.date,'linear','extrap');


nccreate(ficout,'TIME', 'Dimensions',...
    {'TIME' Nhr},'Datatype','single','FillValue',single(99999));
ncwriteatt(ficout,'TIME','units','days since 1950-01-01 00:00:00 UTC');
nccreate(ficout,'LATITUDE', 'Dimensions',...
    {'TIME' },'Datatype','single','FillValue',single(99999));
ncwriteatt(ficout,'LATITUDE','units','degree_north');
nccreate(ficout,'LONGITUDE', 'Dimensions',...
    {'TIME' },'Datatype','single','FillValue',single(99999));
ncwriteatt(ficout,'LONGITUDE','units','degree_east');
nccreate(ficout,'PRES', 'Dimensions',...
    {'TIME' },'Datatype','single','FillValue',single(99999));
ncwriteatt(ficout,'PRES','units','SEA PRESSURE');
nccreate(ficout,'TEMP', 'Dimensions',...
    {'TIME'},'Datatype','single','FillValue',single(99999));
ncwriteatt(ficout,'TEMP','units','degree_Celsius');
nccreate(ficout,'PSAL', 'Dimensions',...
    {'TIME'},'Datatype','single','FillValue',single(99999));
ncwriteatt(ficout,'PSAL','units','PRACTICAL SALINITY');
ncwrite(ficout,'TIME',single(hrdata.date-datenum(1950,1,1)));
ncwrite(ficout,'TEMP',single(hrdata.T));
ncwrite(ficout,'PSAL',single(hrdata.S));
ncwrite(ficout,'PRES',single(hrdata.P));
if hrdata.isfluo
    nccreate(ficout,'CHLA', 'Dimensions',...
        {'TIME'},'Datatype','single','FillValue',single(99999));
    ncwriteatt(ficout,'CHLA','units','mg/m3');
    ncwrite(ficout,'CHLA',single(hrdata.F));
end
if hrdata.isoxy
    nccreate(ficout,'DOXY', 'Dimensions',...
        {'TIME'},'Datatype','single','FillValue',single(99999));
    ncwriteatt(ficout,'DOXY','units','micromole/kg');
    ncwrite(ficout,'DOXY',single(hrdata.O));
end
ncwrite(ficout,'LATITUDE',single(latinterp));
ncwrite(ficout,'LONGITUDE',single(loninterp));


for ii=1:length(listatt),
    k=strfind(listatt{ii},'number');
    if length(k)==0
        ncwriteatt(ficout,'/',listatt{ii},getfield(data_att,listatt{ii}));
    end
end

%
ncwriteatt(ficout,'/','DATE_CREATION',datestr(now(),'yyyymmddHHMMSS'));

name_prof = ficout;
sc_write_global_attribute;


