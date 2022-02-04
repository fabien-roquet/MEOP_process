plot_diags=0;
N1=0; N2=0;

% load data and init flags
Mqc=ARGO_load_qc(name_prof,0);
Mqc.PRES_QC(Mqc.PRES_QC==0) = 1;
Mqc.TEMP_QC(Mqc.TEMP_QC==0) = 1;
Mqc.PSAL_QC(Mqc.PSAL_QC==0) = 1;
if isfield(Mqc,'CHLA_QC'), Mqc.CHLA_QC(Mqc.CHLA_QC==0) = 1; end
if isfield(Mqc,'DOXY_QC'), Mqc.DOXY_QC(Mqc.DOXY_QC==0) = 1; end
if isfield(Mqc,'LIGHT_QC'), Mqc.LIGHT_QC(Mqc.LIGHT_QC==0) = 1; end
ARGO_save_qc(name_prof,Mqc,0);


% create default paramaters
if ~any(strcmp(conf.table_param.Properties.RowNames,EXP)),
    temp_error=0.1; psal_error=0.2; minT=-3; maxT=32; minS= 4; maxS=40; min_Nprof= 30;
    pmax = 1000; pmax_fluo = 200; is_lon_centre_180 = 0;
    conf.table_param(EXP,:)={temp_error psal_error minT maxT minS maxS min_Nprof pmax pmax_fluo is_lon_centre_180};
    name_file=[conf.processdir 'table_param.csv'];
    writetable(conf.table_param,name_file,'WriteRowNames',1,'Delimiter',',');    
end

% max/min paramters
minT=conf.table_param{EXP,'minT'};
maxT=conf.table_param{EXP,'maxT'};
minS=conf.table_param{EXP,'minS'};
maxS=conf.table_param{EXP,'maxS'};
min_Nprof=conf.table_param{EXP,'min_Nprof'};

%% lat/lon/date
[Mqc,nn]=remove_profiles(info_deployment,smru_name,'index',find(isnan(Mqc.LATITUDE.*Mqc.LONGITUDE.*Mqc.JULD)),suffix);
N1=nn;

%% outliers
nT=nansum(double(Mqc.TEMP_QC<=1));
[Mqc,nn1]=remove_profiles(info_deployment,smru_name,'index',find(nT<3&nT>0),suffix);
[Mqc,nn2]=remove_profiles(info_deployment,smru_name,'Tmin',minT,suffix);
[Mqc,nn3]=remove_profiles(info_deployment,smru_name,'Tmax',maxT,suffix);
[Mqc,nn4]=remove_profiles(info_deployment,smru_name,'index',find(all(diff(Mqc.TEMP)==0)),suffix);
[Mqc,nn5]=remove_profiles(info_deployment,smru_name,'index',find(Mqc.LATITUDE==0),suffix);
[Mqc,nn6]=remove_profiles(info_deployment,smru_name,'index',find(abs(Mqc.LATITUDE-nanmean(Mqc.LATITUDE))>5*nanstd(Mqc.LATITUDE)),suffix);
N1=nn1+nn2+nn3+nn4+nn5+nn6;

if nn5, disp(sprintf('Bad locations (lat=0): %d profiles',nn5)), end

%% salinity number point
nS=nansum(double(Mqc.PSAL_QC<2));

[Mqc,nn1]=remove_Sprofiles(info_deployment,smru_name,'index',find(nS<=5&nS>0),suffix);
%[Mqc,nn1]=remove_Sprofiles(info_deployment,smru_name,'index',find(nS<=8&nS>0),suffix);
[Mqc,nn2]=remove_Sprofiles(info_deployment,smru_name,'Smin',minS,suffix);
[Mqc,nn3]=remove_Sprofiles(info_deployment,smru_name,'Smax',maxS,suffix);
%I=find(Mqc.PSAL_QC<=1 & Mqc.TEMP_QC>1); remove_Sdata;
N2=nn1+nn2+nn3;

%% density calculation
nS=nansum(double(Mqc.PSAL_QC<2));
[maxP,ii]=max(Mqc.PRES); Ibot=sub2ind([Mqc.nr Mqc.np],max(ii-1,1),1:Mqc.np);

% vertical resolution
I=find(any(diff(Mqc.PRES)>ones(Mqc.nr-1,1)*maxP/3));
[Mqc,nn1]=remove_Sprofiles(info_deployment,smru_name,'index',I,suffix);

% top-to-bottom inversion
Mqc.PTMP = sw_ptmp( Mqc.PSAL, Mqc.TEMP, Mqc.PRES, 0 );
Mqc.SIG0 = sw_dens0( Mqc.PSAL, Mqc.PTMP ) - 1000 ;
I=find(Mqc.PTMP(2,:)<7 & nS>5 & Mqc.SIG0(Ibot)-Mqc.SIG0(2,:)<-.03);
[Mqc,nn2]=remove_Sprofiles(info_deployment,smru_name,'index',I,suffix);

N2=N2+nn1+nn2;
% add default value for new tag in table_coeff
if sum(strcmp(smru_name,conf.table_coeff.smru_platform_code))==0
    new_tag = {0,0,0,0,0,0,'no comment'};
    conf.table_coeff=[conf.table_coeff;new_tag];
    conf.table_coeff.Properties.RowNames{end}=smru_name;
end
% manual editing
if conf.table_coeff{smru_name,'remove'},
    Mqc=remove_tag(info_deployment,smru_name);
end

if conf.table_coeff{smru_name,'Sremove'},
    [Mqc,nn1]=remove_Sprofiles(info_deployment,smru_name);
end
N2=N2+nn1;

filters1=[];
if any(strcmp(EXP, conf.table_filter.smru_platform_name)),
    filters1 = conf.table_filter(strcmp(EXP, conf.table_filter.smru_platform_name),:);
end

filters2=[];
if any(strcmp(smru_name, conf.table_filter.smru_platform_name)),
    filters2 = conf.table_filter(strcmp(smru_name, conf.table_filter.smru_platform_name),:);
end

filters = [filters1; filters2];
if ~isempty(filters),
    for kk=1:length(filters.smru_platform_name),
        if filters{kk,'Sonly'},
            [Mqc,nn1]=remove_Sprofiles(info_deployment,smru_name,...
                filters{kk,'filter'}{1},filters{kk,{'x1','x2'}},suffix);
            N2=N2+nn1;
        else
            [Mqc,nn1]=remove_profiles(info_deployment,smru_name,...
                filters{kk,'filter'}{1},filters{kk,{'x1','x2'}},suffix);
            N1=N1+nn1;
        end
    end
end
    


% %% spike test 1
% Mqc.PTMP = sw_ptmp( Mqc.PSAL, Mqc.TEMP, Mqc.PRES, 0 );
% Mqc.SIG0 = sw_dens0( Mqc.PSAL, Mqc.PTMP ) - 1000 ;
% 
% if Mqc.np>2,
%     sD2=[zeros(1,Mqc.np)*NaN;abs(diff(Mqc.SIG0,2))/2;zeros(1,Mqc.np)*NaN];
%     dD2=[zeros(1,Mqc.np)*NaN;abs(Mqc.SIG0(3:end,:)-Mqc.SIG0(1:end-2,:))/2;zeros(1,Mqc.np)*NaN ];
%     I=find(sD2>max(4*dD2,.1)); remove_Sdata;
% end


% %% density inversion 2
% Mqc.PTMP = sw_ptmp( Mqc.PSAL, Mqc.TEMP, Mqc.PRES, 0 );
% Mqc.SIG0 = sw_dens0( Mqc.PSAL, Mqc.PTMP ) - 1000 ;
% 
% nn1=0; Iprof=[];
% if plot_diags, figure(10),clf,hold on, end
% for kk=1:Mqc.np,
%     Pint=(0:1000)';
%     if all(isnan(Mqc.SIG0(:,kk))), continue, end
%     ii = find(~isnan(Mqc.SIG0(:,kk))) ;
%     if length(ii)>1
%         Dint=interp1(Mqc.PRES(ii,kk),Mqc.SIG0(ii,kk),Pint);
%         if any(Dint(51:end)-Dint(1:end-50)<-.2)
%             if plot_diags, plot(Mqc.SIG0(:,kk),-Mqc.PRES(:,kk)); end
%             Iprof(end+1)=kk;
%             nn1=nn1+1;
%         end
%     end
% end
% [Mqc,nn2]=remove_Sprofiles(info_deployment,smru_name,'index',Iprof,suffix);
% N2=N2+nn1;


%% minimum number of "min_Nprof [default:30]" valid T profiles

nT=nansum(double(Mqc.TEMP_QC(:,:)<=1));
if length(find(nT>5))<min_Nprof & length(find(nT>5))>0,
    Mqc = remove_tag(info_deployment,smru_name);
    conf.table_coeff{smru_name,'remove'} = 1;
    name_file=[conf.processdir 'table_coeff.csv'];
    writetable(conf.table_coeff,name_file,'WriteRowNames',1,'Delimiter',',');
else
    disp(sprintf('  %s: %d profiles and %d Sprofiles removed',smru_name,N1,N2));    
end


