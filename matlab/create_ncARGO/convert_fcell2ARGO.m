function convert_fcell2ARGO(conf,EXP,name_fcell,suffix,Nlevels,one_smru_name)
% convert_fcell(ficin,ficout)
%   conversion du format fcell au format argo

load(name_fcell);

Nprof=size(hi,1);
Nparam=3+isoxy+isfluo;
if exist('islight','var'),
    Nparam = Nparam + islight;
else
    islight=0;
end

fixed_levels=1;
if ~exist('Nlevels','var') || isempty(Nlevels),
    fixed_levels=0; 
    Nlevels=max(hi(:,8)); 
end

if ~exist('one_smru_name','var')
    one_smru_name = '';
end

date_ref=datenum('19500101000000','yyyymmddHHMMSS');

Juld=zeros(Nprof,1)*NaN;
Lat=zeros(Nprof,1)*NaN;
Lon=zeros(Nprof,1)*NaN;

Pres=zeros(Nlevels,Nprof)*NaN;
Temp=zeros(Nlevels,Nprof)*NaN;
Sali=zeros(Nlevels,Nprof)*NaN;
Fluo=zeros(Nlevels,Nprof)*NaN;
Oxy=zeros(Nlevels,Nprof)*NaN;
Light=zeros(Nlevels,Nprof)*NaN;

Pqc=repmat('9',Nlevels,Nprof);
Tqc=repmat('9',Nlevels,Nprof);
Sqc=repmat('9',Nlevels,Nprof);
Fqc=repmat('9',Nlevels,Nprof);
Oqc=repmat('9',Nlevels,Nprof);
Lqc=repmat('9',Nlevels,Nprof);

platform_number=repmat(' ',8,Nprof);
cycle_number=zeros(Nprof,1)*NaN;
dc_ref=repmat(' ',32,Nprof);

ltag=unique(hi(:,10));n=0;
Ntag=length(ltag);

for tag=1:Ntag,
    
    Itag=find(hi(:,10)==ltag(tag));
    N=length(Itag);
    Juld(1+n:N+n)=hi(Itag,4)-date_ref;
    Lat(1+n:N+n)=hi(Itag,5);
    Lon(1+n:N+n)=hi(Itag,6);
    cycle_number(1+n:N+n)=1:N;
    for ii=1:N,
        kk=Itag(ii);
        if any(ltag)>99, error('too much tags in the dataset (max ntag=100)'); end
        platform_number(:,ii+n)=sprintf('%6s%02d',EXP(1:min(6,length(EXP))),ltag(tag));
        dc_ref(:,ii+n)=sprintf('%24s%02d%06d',EXP(1:min(6,length(EXP))),ltag(tag),ii);
        Np=size(PTi{kk},1);
        Ns=size(PSi{kk},1);
        if Np>Nlevels,
            disp(['N>Nlevels ' num2str(Np)]);
        end
        if Np>0
            Pres(1:Np,ii+n)=PTi{kk}(:,1)'; Pqc(                1:Np,ii+n)='0';
            Temp(1:Np,ii+n)=PTi{kk}(:,2)'; Tqc(~isnan(PTi{kk}(:,2)),ii+n)='0';
        end
        if Ns>0
            Sali(1:Ns,ii+n)=PSi{kk}(:,2)'; Sqc(~isnan(PSi{kk}(:,2)),ii+n)='0';
        end
        if isfluo & Np>0
            Fluo(1:Np,ii+n)=PFi{kk}(:,2)'; Fqc(~isnan(PFi{kk}(:,2)),ii+n)='0';
        end
        if isoxy & Np>0
            Oxy (1:Np,ii+n)=POi{kk}(:,2)'; Oqc(~isnan(POi{kk}(:,2)),ii+n)='0';
        end
        if islight & Np>0
            Light (1:Np,ii+n)=PLi{kk}(:,2)'; Lqc(~isnan(PLi{kk}(:,2)),ii+n)='0';
        end
    end
    n=n+N;
    
end

%% create netcdf-argo

platform_number_cell=cellstr(platform_number(1:8,:)');
list_tag = unique(platform_number_cell);

for ii=1:length(list_tag),
    
    I=find(strcmp(list_tag{ii},platform_number_cell));
    if ~fixed_levels,
        Nlevels=max([sum(any(~isnan(Temp(:,I)),2)) sum(any(~isnan(Sali(:,I)),2))]);
    end
    if ~Nlevels,
        disp(['No valid TS data'])
        continue,
    end
    
    smru_name = char(hs(I(1)));
    [smru_prefix,Nsplit] = Nsplit_from_smru_name(smru_name);
    [one_smru_prefix,Nsplit] = Nsplit_from_smru_name(one_smru_name);
    if ~isempty(one_smru_name) & ~strcmp(one_smru_prefix,smru_prefix),
        continue
    end
    K=find(strcmp(conf.list_smru_platform_code,smru_prefix));
    
    if length(K)>0
        
        Nprof=length(I);
        ficoutind = sprintf('%s%s/%s_%s',conf.datadir,EXP,smru_name,suffix);
        
        if exist(ficoutind), delete(ficoutind); end
        if ~exist('isfluo','var'), isfluo = 0; end
        if ~exist('isoxy','var'), isoxy = 0; end
        if ~exist('islight','var'), islight = 0; end
        
        ARGO_create(ficoutind,Nprof,Nlevels,isfluo,isoxy,islight,0);
        %ncwrite(ficoutind,'PLATFORM_NUMBER',platform_number(:,I));
        ncwrite(ficoutind,'PLATFORM_NUMBER',...
            repmat(sprintf('%08d',str2num(conf.platform_json(K).platform_code)),Nprof,1)');
        ncwrite(ficoutind,'PI_NAME',repmat(sprintf('%64s',PI),Nprof,1)');
        ncwrite(ficoutind,'PROJECT_NAME',repmat(sprintf('%64s','MEOP'),Nprof,1)');
        ncwrite(ficoutind,'CYCLE_NUMBER',cycle_number(I));
        
        ncwrite(ficoutind,'JULD',Juld(I));
        ncwrite(ficoutind,'JULD_LOCATION',Juld(I));
        ncwrite(ficoutind,'LATITUDE', Lat(I));
        ncwrite(ficoutind,'LONGITUDE', Lon(I));
        ncwrite(ficoutind,'PRES', Pres(1:Nlevels,I));
        ncwrite(ficoutind,'PRES_QC',Pqc(1:Nlevels,I));
        ncwrite(ficoutind,'TEMP', Temp(1:Nlevels,I));
        ncwrite(ficoutind,'TEMP_QC', Tqc(1:Nlevels,I));
        ncwrite(ficoutind,'PSAL', Sali(1:Nlevels,I));
        ncwrite(ficoutind,'PSAL_QC', Sqc(1:Nlevels,I));
        if isfluo,
            ncwrite(ficoutind,'CHLA', Fluo(1:Nlevels,I));
            ncwrite(ficoutind,'CHLA_QC', Fqc(1:Nlevels,I));
        end
        if isoxy,
            ncwrite(ficoutind,'DOXY', single(Oxy(1:Nlevels,I)));
            ncwrite(ficoutind,'DOXY_QC', Oqc(1:Nlevels,I));
        end
        if islight,
            ncwrite(ficoutind,'LIGHT', single(Light(1:Nlevels,I)));
            ncwrite(ficoutind,'LIGHT_QC', Lqc(1:Nlevels,I));
        end
        ncwriteatt(ficoutind,'/','comment',' ');
        
    end
    
end




