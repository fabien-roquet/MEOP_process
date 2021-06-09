function ARGO_save_qc(name_prof,argo_qc,adjusted);

if ~adjusted,
    
    a=argo_qc.PRES_QC;
    for kk=0:9,
        a(a==kk)=48+kk;
    end
    data = char(a);
    ncwrite(name_prof,'PRES_QC',data);
    ncwrite(name_prof,'PRES_ADJUSTED_QC',data);
    
    a=argo_qc.TEMP_QC;
    for kk=0:9,
        a(a==kk)=48+kk;
    end
    data = char(a);
    ncwrite(name_prof,'TEMP_QC',data);
    ncwrite(name_prof,'TEMP_ADJUSTED_QC',data);
    
    a=argo_qc.PSAL_QC;
    for kk=0:9,
        a(a==kk)=48+kk;
    end
    data = char(a);
    ncwrite(name_prof,'PSAL_QC',data);
    ncwrite(name_prof,'PSAL_ADJUSTED_QC',data);
    
    if isfield(argo_qc,'CHLA')
        a=argo_qc.CHLA_QC;
        for kk=0:9,
            a(a==kk)=48+kk;
        end
        data = char(a);
        try
            ncwrite(name_prof,'CHLA_QC',data);
            ncwrite(name_prof,'CHLA_ADJUSTED_QC',data);
        end
    end
    
    if isfield(argo_qc,'DOXY')
        a=argo_qc.DOXY_QC;
        for kk=0:9,
            a(a==kk)=48+kk;
        end
        data = char(a);
        try
            ncwrite(name_prof,'DOXY_QC',data);
            ncwrite(name_prof,'DOXY_ADJUSTED_QC',data);
        end
    end
    
     if isfield(argo_qc,'LIGHT')
        a=argo_qc.LIGHT_QC;
        for kk=0:9,
            a(a==kk)=48+kk;
        end
        data = char(a);
        try
            ncwrite(name_prof,'LIGHT_QC',data);
            ncwrite(name_prof,'LIGHT_ADJUSTED_QC',data);
        end
    end
    
else
    
    a=argo_qc.PRES_QC;
    for kk=0:9,
        a(a==kk)=48+kk;
    end
    data = char(a);
    ncwrite(name_prof,'PRES_ADJUSTED_QC',data);
    
    a=argo_qc.TEMP_QC;
    for kk=0:9,
        a(a==kk)=48+kk;
    end
    data = char(a);
    ncwrite(name_prof,'TEMP_ADJUSTED_QC',data);
    
    a=argo_qc.PSAL_QC;
    for kk=0:9,
        a(a==kk)=48+kk;
    end
    data = char(a);
    ncwrite(name_prof,'PSAL_ADJUSTED_QC',data);
    
    if isfield(argo_qc,'CHLA')
        a=argo_qc.CHLA_QC;
        for kk=0:9,
            a(a==kk)=48+kk;
        end
        data = char(a);
        try, ncwrite(name_prof,'CHLA_ADJUSTED_QC',data); end
    end
    
    if isfield(argo_qc,'DOXY')
        a=argo_qc.DOXY_QC;
        for kk=0:9,
            a(a==kk)=48+kk;
        end
        data = char(a);
        try, ncwrite(name_prof,'DOXY_ADJUSTED_QC',data); end
    end
    
    if isfield(argo_qc,'LIGHT')
        a=argo_qc.LIGHT_QC;
        for kk=0:9,
            a(a==kk)=48+kk;
        end
        data = char(a);
        try, ncwrite(name_prof,'LIGHT_ADJUSTED_QC',data); end
    end
    
end

%% update descriptive attributes

NS=length(find(nansum((argo_qc.PSAL_QC==1).*(argo_qc.TEMP_QC==1))~=0));
NT=length(find(nansum(argo_qc.TEMP_QC==1)~=0));
ncwriteatt(name_prof,'/','number_of_ts_profiles',NS);
ncwriteatt(name_prof,'/','number_of_t_profiles',NT);

NF=0;
if isfield(argo_qc,'CHLA')
    NF=length(find(nansum(argo_qc.CHLA_QC==1)~=0));
end
ncwriteatt(name_prof,'/','number_chla_profiles',NF);

NO=0;
if isfield(argo_qc,'DOXY')
    NO=length(find(nansum(argo_qc.DOXY_QC==1)~=0));
end
ncwriteatt(name_prof,'/','number_doxy_profiles',NO);

NO=0;
if isfield(argo_qc,'LIGHT')
    NO=length(find(nansum(argo_qc.LIGHT_QC==1)~=0));
end
ncwriteatt(name_prof,'/','number_light_profiles',NO);

Lat=argo_qc.LATITUDE;
Lon=argo_qc.LONGITUDE;
I=find(nansum((argo_qc.PSAL_QC==1)+(argo_qc.TEMP_QC==1))~=0);

if isempty(I),

    ncwriteatt(name_prof,'/','geospatial_lat_min',' ');
    ncwriteatt(name_prof,'/','geospatial_lat_max',' ');
    ncwriteatt(name_prof,'/','geospatial_lon_min',' ');
    ncwriteatt(name_prof,'/','geospatial_lon_max',' ');

else
    
    ncwriteatt(name_prof,'/','geospatial_lat_min',min(Lat(I)));
    ncwriteatt(name_prof,'/','geospatial_lat_max',max(Lat(I)));
    
    mLon=min(Lon(I)); MLon=max(Lon(I));
    mLon(mLon<0)=mLon(mLon<0)+360;
    MLon(MLon<0)=MLon(MLon<0)+360;
    ncwriteatt(name_prof,'/','geospatial_lon_min',mLon);
    ncwriteatt(name_prof,'/','geospatial_lon_max',MLon);
end
