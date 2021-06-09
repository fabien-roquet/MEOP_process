function Mtemp=extract_profil(Mgroup,I)
%% Extract profil matrice type Mgroup
% Mgroup = matrice with one or several data tags plus WOD associated
% I = profil number selected
Mtemp=Mgroup;
if isfield(Mgroup,'LATITUDE')
    Mtemp.LATITUDE=Mtemp.LATITUDE(I);
end
if isfield(Mgroup,'LONGITUDE')
    Mtemp.LONGITUDE=Mtemp.LONGITUDE(I);
end
if isfield(Mgroup,'JULD_LOCATION')
    Mtemp.JULD_LOCATION=Mtemp.JULD_LOCATION(I);
end
if isfield(Mgroup,'CYCLE_NUMBER')
    Mtemp.CYCLE_NUMBER=Mtemp.CYCLE_NUMBER(I);
end
if isfield(Mgroup,'PLATFORM_NUMBER')
    Mtemp.PLATFORM_NUMBER=Mtemp.PLATFORM_NUMBER(:,I);
end
if isfield(Mgroup,'PI_NAME')
    Mtemp.PI_NAME=Mtemp.PI_NAME(:,I);
end
if isfield(Mgroup,'JULD')
    Mtemp.JULD=Mtemp.JULD(I);
end
if isfield(Mgroup,'PRES')
    Mtemp.PRES=Mtemp.PRES(:,I);
end
if isfield(Mgroup,'TEMP')
    Mtemp.TEMP=Mtemp.TEMP(:,I);
end
if isfield(Mgroup,'PSAL')
    Mtemp.PSAL=Mtemp.PSAL(:,I);
end
if isfield(Mgroup,'CHLA')
    Mtemp.CHLA=Mtemp.CHLA(:,I);
end
if isfield(Mgroup,'DOXY')
    Mtemp.DOXY=Mtemp.DOXY(:,I);
end
if isfield(Mgroup,'PRES_QC')
    Mtemp.PRES_QC=Mtemp.PRES_QC(:,I);
end
if isfield(Mgroup,'TEMP_QC')
    Mtemp.TEMP_QC=Mtemp.TEMP_QC(:,I);
end
if isfield(Mgroup,'PSAL_QC')
    Mtemp.PSAL_QC=Mtemp.PSAL_QC(:,I);
end
if isfield(Mgroup,'CHLA_QC')
    Mtemp.CHLA_QC=Mtemp.CHLA_QC(:,I);
end
if isfield(Mgroup,'DOXY_QC')
    Mtemp.DOXY_QC=Mtemp.DOXY_QC(:,I);
end
if isfield(Mgroup,'np')
    Mtemp.np=length(I);
end
if isfield(Mgroup,'platform_number')
    Mtemp.platform_number=Mtemp.platform_number(I);
end
if isfield(Mgroup,'ntag')
    Mtemp.ntag=1;
end
if isfield(Mgroup,'nprof')
    Mtemp.nprof=length(I);
end
if isfield(Mgroup,'index_tag')
    Mtemp.index_tag=Mtemp.index_tag(I);
end
if isfield(Mgroup,'index_tag')
    Mtemp.index_tag(:)=1;
end
if isfield(Mgroup,'platform_number')
    try
    F=find(strcmp(Mgroup.platform_number(I(1)),Mtemp.list_platform_number));
    Mtemp.list_descr=Mtemp.list_platform_number(F);
    Mtemp.list_platform_number = Mtemp.list_platform_number(F);
    Mtemp.list_smru_platform = Mtemp.list_smru_platform(F);
    end
end
if isfield(Mgroup,'geospatial_lat_min')
    Mtemp.geospatial_lat_min= Mtemp.geospatial_lat_min;
end
if isfield(Mgroup,'geospatial_lat_max')
    Mtemp.geospatial_lat_max= Mtemp.geospatial_lat_max;
end
if isfield(Mgroup,'geospatial_lon_min')
    Mtemp.geospatial_lon_min= Mtemp.geospatial_lon_min;
end
if isfield(Mgroup,'geospatial_lon_max')
    Mtemp.geospatial_lon_max= Mtemp.geospatial_lon_max;
end
if isfield(Mgroup,'geospatial_vertical_min')
    Mtemp.geospatial_vertical_min= Mtemp.geospatial_vertical_min;
end
if isfield(Mgroup,'number_of_ts_profiles')
    Mtemp.number_of_ts_profiles=Mtemp.number_of_ts_profiles;
end
if isfield(Mgroup,'number_of_t_profiles')
    Mtemp.number_of_t_profiles=Mtemp.number_of_t_profiles;
end
if isfield(Mgroup,'Pmask')
    Mtemp.Pmask=Mtemp.Pmask(:,I);
end
if isfield(Mgroup,'Tmask')
    Mtemp.Tmask=Mtemp.Tmask(:,I);
end
if isfield(Mgroup,'Smask')
    Mtemp.Smask=Mtemp.Smask(:,I);
end
