function data = load_WOD(conf)
% load WOD data in Argo netCDF format
%
% conf is a struct variable with the following fields
%    woddir : directory of WOD files
%    lat and lon : position of data needed. A WOD data tile is loaded if at least
%        one of the position is inside its area.


data=[];
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
                name_meop=sprintf('%sMEOP-CTD_lon%02d_lat%02d_prof.nc',conf.meopdir,ii,jj);
                if exist(name_meop,'file')
                    aux=ARGO_load(name_meop);
                    if isempty(data), data = aux;
                    else; data = ARGO_concat(data,aux);
                    end
                end
            end
        end
        
    end
    
end
if ~isempty(data)
    data.LONGITUDE(data.LONGITUDE>180)=data.LONGITUDE(data.LONGITUDE>180)-360;
end


