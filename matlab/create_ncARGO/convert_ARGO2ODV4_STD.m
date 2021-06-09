function convert_ARGO2ODV4(argo_data,path)

warning off

PRES=argo_data.PRES;
TEMP=argo_data.TEMP;
SALI=argo_data.PSAL;
Pmask=double(argo_data.PRES_QC==1);
Tmask=double(argo_data.TEMP_QC==1);
Smask=double(argo_data.PSAL_QC==1);
index=Tmask+Smask;
if isempty(index) || sum(index(:))==0, return, end

PRES(Pmask==0) = NaN;
TEMP(Tmask==0) = NaN;
SALI(Smask==0) = NaN;

std_lev = [0:5:95 100:25:475 500:50:2000]'; 
sc_profiles_interp;
Pmask=double(~isnan(PRES_std_lev));
Tmask=double(~isnan(TEMP_std_lev));
Smask=double(~isnan(SALI_std_lev));
index=Tmask+Smask;

cruise=argo_data.smru_platform_code;
I=1:length(argo_data.JULD);%find(strcmp(argo_data.smru_name,cruise));

fid=fopen(path,'w');
fprintf(fid,'// created: %s\n',datestr(now));
fprintf(fid,'%s\n','Cruise	Station	Type	mon/day/yr	hh:mm	Longitude [degrees_east]	Latitude [degrees_north]	Depth [m]	QF	Temperature [°C]	QF	Salinity [psu]	QF');
for ii=1:length(I),
    if sum(Tmask(:,I(ii)))==0, continue, end
    station=sprintf('%04d',ii);
    date=datestr(argo_data.JULD(I(ii)),'mm/dd/yyyy\tHH:MM');
    lon=argo_data.LONGITUDE(I(ii));    lat=argo_data.LATITUDE(I(ii));
    for pp=1:N_std,
        if index(pp,I(ii))==0, continue, end;
        fprintf(fid,'%s\t%s\tC\t%s\t%6.3f\t%6.3f\t%4.1f\t0\t',cruise,station,date,lon,lat,PRES_std_lev(pp,I(ii)));
        if Tmask(pp,I(ii)),
            fprintf(fid,'%8.4f\t0\t',TEMP_std_lev(pp,I(ii)));
        else
            fprintf(fid,'\t1\t');
        end
        if Smask(pp,I(ii)),
            fprintf(fid,'%8.4f\t0\t',SALI_std_lev(pp,I(ii)));
        else
            fprintf(fid,'\t1\t');
        end
        fprintf(fid,'\n');
    end
end
fclose(fid);

