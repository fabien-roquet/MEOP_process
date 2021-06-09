function argo = load_CORA(conf)
% load CORA data in Argo netCDF format
%
% conf is a struct variable with the following fields
%    woddir : directory of WOD files
%    lat and lon : position of data needed. A WOD data tile is loaded if at least
%        one of the position is inside its area.


argo=[];
if isempty(conf),
    return
end

conf.lon(conf.lon>180)=conf.lon(conf.lon>180)-360;
if ~exist(conf.coradir,'dir')
    error('Unable to find location of CORA files')
end

for lat_min=-80:10:70,
    
    lat_limit = [lat_min lat_min+10];
    if any(conf.lat>lat_limit(1) & conf.lat<=lat_limit(2))
        
        for lon_min=-180:10:170,
            lon_limit = [lon_min lon_min+10];
            if any(conf.lon>lon_limit(1) & conf.lon<=lon_limit(2))
                strlon = 'E'; if lon_min<0, strlon = 'W'; end
                strlat = 'N'; if lat_min<0, strlat = 'S'; end
                name_file = sprintf('%sCORA_lon%02d%s_lat%02d%s.nc',conf.coradir,...
                    abs(lon_min),strlon,abs(lat_min),strlat);
                if exist(name_file,'file') && ~isempty(getfield(ncinfo(name_file),'Dimensions'))
                    aux=ARGO_load(name_file);
                    if isempty(argo), argo = aux;
                    else; argo = ARGO_concat(argo,aux);
                    end
                end
            end
        end
        
    end
    
end
if ~isempty(argo)
    argo.LONGITUDE(argo.LONGITUDE>180)=argo.LONGITUDE(argo.LONGITUDE>180)-360;
end

% if any(conf.lat>-43 & conf.lat<-30 & conf.lon>120 & conf.lat<155)
%     load([conf.woddir '../CTD_imos_profil.mat']);
%     argo=ARGO_concat(argo,data);
%     load([conf.woddir '../Argos_imos_profil.mat']);
%     argo=ARGO_concat(argo,data);
% end

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

