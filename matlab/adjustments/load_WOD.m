function argo_wod = load_WOD(conf)
% load WOD data in Argo netCDF format
%
% conf is a struct variable with the following fields
%    woddir : directory of WOD files
%    lat and lon : position of data needed. A WOD data tile is loaded if at least
%        one of the position is inside its area.

argo_wod=[];
if isempty(conf),
    return
end

conf.lon(conf.lon>180)=conf.lon(conf.lon>180)-360;

for jj=1:8,
    
    lat_limit = [0 20]+(jj-5)*20;
    if any(conf.lat>lat_limit(1) & conf.lat<lat_limit(2))
        
        for ii=1:12,
            lon_limit = [0 30]+(ii-1)*30;
            if lon_limit(1)>181, lon_limit=lon_limit-360; end
            if any(conf.lon>lon_limit(1) & conf.lon<lon_limit(2))
                name_wod=sprintf('%sWOD_lon%02d_lat%02d_prof.nc',conf.woddir,ii,jj);
                if exist(name_wod,'file')
                    aux=ARGO_load(name_wod);
                    if isempty(argo_wod), argo_wod = aux;
                    else; argo_wod = ARGO_concat(argo_wod,aux);
                    end
                end
            end
        end
        
    end
    
end

if ~isempty(argo_wod)
    argo_wod.LONGITUDE(argo_wod.LONGITUDE>180)=argo_wod.LONGITUDE(argo_wod.LONGITUDE>180)-360;
end

if any(conf.lat>-43 & conf.lat<-30 & conf.lon>120 & conf.lat<155)
    load([conf.woddir '../CTD_imos_profil.mat']);
    argo_wod=ARGO_concat(argo_wod,data);
    load([conf.woddir '../Argos_imos_profil.mat']);
    argo_wod=ARGO_concat(argo_wod,data);
end

% if (min(conf_clim.lat)>=-50 & max(conf_clim.lat)<=-30 & min(conf_clim.lon)>=120 & max(conf_clim.lon)<=50)
%     dir_ctd_imos='../CALIBRATION/SARDI_data/IMOS_CTD/IMOS_-_Australian_National_Mooring_Network_(ANMN)_-_CTD_Profiles/';
%     list_file=dir([dir_ctd_imos '*.nc']);
%     data=[];
%     for ii=1:length(list_file)
%         aux=[];
%         aux.PRES=ncread([dir_ctd_imos list_file(ii).name],'DEPTH');
%         aux.LATITUDE=ncread([dir_ctd_imos list_file(ii).name],'LATITUDE');
%         aux.LONGITUDE=ncread([dir_ctd_imos list_file(ii).name],'LONGITUDE');
%         aux.JULD=ncread([dir_ctd_imos list_file(ii).name],'TIME')+datenum(1950,1,1);
%         aux.TEMP=ncread([dir_ctd_imos list_file(ii).name],'TEMP');
%         aux.PSAL=ncread([dir_ctd_imos list_file(ii).name],'PSAL');
%         aux.platform_number{1}=ncreadatt([dir_ctd_imos list_file(ii).name],'/','instrument_serial_number');
%         if size(aux.TEMP,1)==1
%             aux.TEMP=aux.TEMP';
%         end
%         if size(aux.PSAL,1)==1
%             aux.PSAL=aux.PSAL';
%         end
%         aux.np=size(aux.PRES,2);
%         aux.nr=size(aux.PRES,1);
%         data = ARGO_concat(data,aux);
%
%     end
% end

