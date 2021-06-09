%%%%%%%%%%%%%%Correct salinity and Temperature using sal_cor.m and Barker - July 2017

clc
clf
%close all;
%% PARAMETERS TO CHOOSE 

%plot tous les profils 1 ou 0
all_plot=1;
% time threshold to match profiles 0.0069 is 10 min %0.014is 20 min, 0.01 is 15 min, 0.021 is 30 min , 0.03125 % corresponds to 45min in julian day 
thres_min = 0.0069;
%%
format loose
format compact


MqHR=ARGO_load_qc(name_prof,1);
% MqLR=ARGO_load_qc(nc_nameLR,1);
nc_att=ncloadatt_struct(name_prof);

%Load High Resolution data
lat  = MqHR.LATITUDE(1:end);
lon  = MqHR.LONGITUDE(1:end);
Date  = MqHR.JULD_LOCATION(1:end);
Temp  = MqHR.TEMP(:,1:end);
Sal  = MqHR.PSAL(:,1:end);
Pres  = MqHR.PRES(:,1:end);

%Load Low Resolution data
% latL  = ncread(nc_nameLR,'LATITUDE');
% lonL  = ncread(nc_nameLR,'LONGITUDE');
% TempL  = ncread(nc_nameLR,'TEMP_ADJUSTED');
% SalL  = ncread(nc_nameLR,'PSAL_ADJUSTED');
% PresL  = ncread(nc_nameLR,'PRES_ADJUSTED');
% DateL  = ncread(nc_nameLR,'JULD_LOCATION');
%% quality check
%Make NaN missing values
% latL(find(latL==99999))  = NaN;
% lonL(find(lonL==99999))  = NaN;
% TempL(find(TempL==99999))  = NaN;
% SalL(find(SalL==99999))  = NaN;
% PresL(find(PresL==99999))  = NaN;
% DateL(find(DateL==99999))  = NaN;

% % remplace bad quality salinity values by NaN 
% psal_qc=ncread(nc_nameHR,'PSAL_ADJUSTED_QC');
% [x y]=find(psal_qc=='4');
% y=unique(y);

%% Match profiles with time threshold  (thres_min)
% ind=1;
% for k = 1:size(TempL,2)
%     [m,DP]=min(abs((Date-DateL(k))));
%     b(k)=m;
%     if m<thres_min
%         d_date(ind)=Date(DP);
%         ind_hr(ind)=DP;
%         ind_lr(ind)=k;
%         ind=ind+1;
%     end
% end

%% removed wrong profile-pairs based on raw data quality check 

% remplace bad quality salinity values by NaN 
% psal_qc=ncread(nc_nameHR,'PSAL_ADJUSTED_QC');
% [x y]=find(psal_qc=='4');
% y=unique(y);
% y=intersect(y,ind_hr);
% ind_rem=[];
% for i=1:size(y,1)
%     ind_rem(i)=find(ind_hr==y(i));
% end
% 
% disp(['Numbers of pairs (LR+HR) of profiles removed : ',num2str(size(y,1)) , '/', num2str(size(ind_lr,2)), ', based on the HR PSAL quality flag (flag equals to 4)'])
% ind_hr(ind_rem)=[];
% ind_lr(ind_rem)=[];
% d_date(ind_rem)=[];

%% apply step 1 to selected pairs of matching profiles: Thermal Mass Effect - Sal_cor with alpha=beta=0.03 (Vigan's coefficients)
Pi=(0:1000)';
if hr==0
    TiL=repmat(Pi,1,size(Temp,2))*NaN;
    SiL=repmat(Pi,1,size(Temp,2))*NaN;
    PiL=repmat(Pi,1,size(Temp,2))*NaN;

    for jj=1:size(Temp,2)
        I = find(~isnan(Pres(:,jj)));
        if length(I)>1
            TiL(:,jj)=interp1(Pres(I,jj),Temp(I,jj),Pi);
            SiL(:,jj)=interp1(Pres(I,jj),Sal(I,jj),Pi);
            PiL(:,jj)=Pi;
        else
            TiL(:,jj)=NaN;
            SiL(:,jj)=NaN;
            PiL(:,jj)=Pi;
            
        end
    end
    Temp=[];
    Temp=TiL;
    Sal=[];
     Sal=SiL;
     Pres=[];
     Pres=PiL;
     clear SiL TiL PiL
end

TempHR_cor=Temp*NaN;
SalHR_cor=Temp*NaN;
DensHR_cor=Temp*NaN;
TempHR_brut=Temp*NaN;
SalHR_brut=Temp*NaN;
DensHR_brut=Temp*NaN;
SAHR_cor=Temp*NaN;
CTHR_cor=Temp*NaN;
SAHR_brut=Temp*NaN;
CTHR_brut=Temp*NaN;


    
    
% TempLR_brut=[];
% TempLR_cor=[];
% SalLR_cor=[];
% SalLR_brut=[];
% DensLR_cor=[];
% DensLR_brut=[];
% SALR_cor=[];
% CTLR_cor=[];
% SALR_brut=[];
% CTLR_brut=[];

nbpt_brut=[];
d_date=[];

alpha=0.03;
beta=0.03;
SI=1;
a=(2*alpha)./((SI*beta)+2);
b=1-((2*a)/alpha);
disp(['TME parameters: SI=1, alpha = beta = 0.03, a = ',num2str(a),' and b = ', num2str(b)])

for k = 1:size(Temp,2)
        %% correction HR
        T = Temp(:,k);
        S = Sal(:,k);
        p = Pres(:,k);
        d_date=Date(k);
         
        saHR_b=gsw_SR_from_SP(S); % Absolute salinity (SA) from practical salinity (SP)
        ctHR_b=gsw_CT_from_t(saHR_b,T,p); % Conservative Temperature from in-situ temperature
        D=gsw_sigma0(saHR_b,ctHR_b);
        
        [Sc Tc] = sal_cor(alpha,beta,SI,T,S,p);
 
        sa=gsw_SR_from_SP(Sc); % Absolute salinity (SA) from practical salinity (SP)
        ct=gsw_CT_from_t(sa,T,p); % Conservative Temperature from in-situ temperature
        Dc=gsw_sigma0(sa,ct);
         
        SAHR_cor(:,k)=sa;
        CTHR_cor(:,k)=ct;
        TempHR_cor(:,k)=Tc;
        SalHR_cor(:,k)=Sc;
        DensHR_cor(:,k)=Dc;
        TempHR_brut(:,k)=T;
        SalHR_brut(:,k)=S;
        DensHR_brut(:,k)=D;
        SAHR_brut(:,k)=saHR_b;
        CTHR_brut(:,k)=ctHR_b;
        
        % si on est en BR on interpolle sur la grille à 1m

        
        %% correction LR
%         TL = TempL(:,ind_lr(k));
%         SL = SalL(:,ind_lr(k));
%         pL = PresL(:,ind_lr(k));
%        % ind_lr(ind)=k;
%         
% %         sa_trans = gsw_SA_from_SP(SL,pL,lonL(k),latL(k));
% %         ct_trans = gsw_CT_from_t(sa_trans,TL,pL);
% %         dens_trans = gsw_sigma0(sa_trans,ct_trans);
% 
%         I = find(~isnan(pL));
%         nbpt_brut(k,1)=length(I);
%         Pi=(0:1000)';
%         TiL=interp1(pL(I),TL(I),Pi);
%         SiL=interp1(pL(I),SL(I),Pi);
% 
%         
%         saLR_b=gsw_SA_from_SP(SiL,Pi,lonL(ind_lr(k)),latL(ind_lr(k))); % Absolute salinity (SA) from practical salinity (SP)
%         %saLR_b = interp1(pL(I),sa_trans(I),Pi);
%         ctLR_b=gsw_CT_from_t(saLR_b,TiL,Pi); % Conservative Temperature from in-situ temperature
%         DL=gsw_sigma0(saLR_b,ctLR_b);
%         %DL = interp1(pL(I),dens_trans(I),Pi);
%         
%         [ScL TcL] = sal_cor_20161214(alpha,beta,SI,TiL,SiL,Pi);
% 
%         saL=gsw_SA_from_SP(ScL,Pi,lonL(ind_lr(k)),latL(ind_lr(k))); % Absolute salinity (SA) from practical salinity (SP)
%         ctL=gsw_CT_from_t(saL,TiL,Pi); % Conservative Temperature from in-situ temperature
%         DcL=gsw_sigma0(saL,ctL);
% 
%         SALR_cor(:,k)=saL;
%         CTLR_cor(:,k)=ctL;
%         TempLR_cor(:,k)=TcL;
%         SalLR_cor(:,k)=ScL;
%         DensLR_cor(:,k)=DcL;
%         SALR_brut(:,k)=saLR_b;
%         CTLR_brut(:,k)=ctLR_b;
%         TempLR_brut(:,k)=TiL;
%         SalLR_brut(:,k)=SiL;
%         DensLR_brut(:,k)=DL;
       
end
clear saHR_b ctHR_b D T S Sc p sa ct Dc saLR_b ctLR_b TL SL ScL SiL pL saL ctL DcL DP k m ind I Tc TcL TiL DL A B count

%% Step #2 - Density Inversion Removal (Paul Barker and Trevor McDougall, 2017)
SAHR_2=gsw_stabilise_SA_const_t(SAHR_cor,CTHR_cor,Pres);
CTHR_2=gsw_CT_from_t(SAHR_2,TempHR_brut,Pres);
DensHR_2=gsw_sigma0(SAHR_2,CTHR_2);
% 
% SALR_2=gsw_stabilise_SA_const_t(SALR_cor,CTLR_cor,Pres,lonL(ind_lr),latL(ind_lr));
% CTLR_2=gsw_CT_from_t(SALR_2,TempLR_brut,Pres(:,ind_lr));
% DensLR_2=gsw_sigma0(SALR_2,CTLR_2);

%% Step #3 - Gaussian filter smoothing 
firstnan_ctHR=find_ndim(~isnan(CTHR_2),1,'first');
lastnan_ctHR=find_ndim(~isnan(CTHR_2),1,'last');

firstnan_saHR=find_ndim(~isnan(SAHR_2),1,'first');
lastnan_saHR=find_ndim(~isnan(SAHR_2),1,'last');


% firstnan_ctLR=find_ndim(~isnan(CTLR_2),1,'first');
% lastnan_ctLR=find_ndim(~isnan(CTLR_2),1,'last');
% 
% firstnan_saLR=find_ndim(~isnan(SALR_2),1,'first');
% lastnan_saLR=find_ndim(~isnan(SALR_2),1,'last');

lastnan_densHR=find_ndim(~isnan(DensHR_2),1,'last');
%% find if a nan value is in between non nan-values and set to nan all values in between NaNs so that the gaussian filter doesnt propagate NaNs more than necessary
for i=1:size(SAHR_2,2)
    if firstnan_saHR(i)~=0 & lastnan_saHR(i)~=0
        if sum(isnan(SAHR_2(firstnan_saHR(i):lastnan_saHR(i),i)))~= 0
            x=find(isnan(SAHR_2(firstnan_saHR(i):lastnan_saHR(i),i)));
            row=max(x);
            SAHR_2(1:x,i)=NaN;
        end
    end
            
end
for i=1:size(CTHR_2,2)
    if firstnan_ctHR(i)~=0 & lastnan_ctHR(i)~=0
        if sum(isnan(CTHR_2(firstnan_ctHR(i):lastnan_ctHR(i),i)))~= 0
            x=find(isnan(CTHR_2(firstnan_ctHR(i):lastnan_ctHR(i),i)));
            row=max(x);
            CTHR_2(1:x,i)=NaN;
        end
    end
end
for i=1:size(DensHR_2,2)
    if firstnan_saHR(i)~=0 & lastnan_saHR(i)~=0
        if sum(isnan(SAHR_2(firstnan_saHR(i):lastnan_saHR(i),i)))~= 0
            x=find(isnan(SAHR_2(firstnan_saHR(i):lastnan_saHR(i),i)));
            row=max(x);
            DensHR_2(1:x,i)=NaN;
        end
    end
end
firstnan_saHR=find_ndim(~isnan(SAHR_2),1,'first');
firstnan_ctHR=find_ndim(~isnan(CTHR_2),1,'first');
firstnan_densHR=find_ndim(~isnan(DensHR_2),1,'first');
%% ordre 1 
SAHR_3=SAHR_2*NaN;
CTHR_3=CTHR_2*NaN;

% SALR_3=SALR_2*NaN;
% CTLR_3=CTLR_2*NaN;

win=1; % 1dbar window
num=win*5;

for i=1:size(SAHR_2,2)
    
    if firstnan_saHR(i) ~= 0
        c=SAHR_2(firstnan_saHR(i):lastnan_saHR(i),i);
        
        a=repmat(c(1),num,1);
        b=repmat(c(end),num,1);
        d=[a;c;b];
        s=gaussfilt(1:length(d),d,win);
        SAHR_3(firstnan_saHR(i):lastnan_saHR(i),i)=s(num+1:end-num);
    end
        
    if firstnan_ctHR(i) ~= 0
        c=CTHR_2(firstnan_ctHR(i):lastnan_ctHR(i),i);
        
        a=repmat(c(1),num,1);
        b=repmat(c(end),num,1);
        d=[a;c;b];
        f=gaussfilt(1:length(d),d,win);
        CTHR_3(firstnan_ctHR(i):lastnan_ctHR(i),i)=f(num+1:end-num);
    end
    
%     if firstnan_saLR(i) ~= 0
%         c=SALR_2(firstnan_saLR(i):lastnan_saLR(i),i);
%         
%         a=repmat(c(1),num,1);
%         b=repmat(c(end),num,1);
%         d=[a;c;b];
%         s=gaussfilt(1:length(d),d,win);
%         SALR_3(firstnan_saLR(i):lastnan_saLR(i),i)=s(num+1:end-num);
%     end
        
%     if firstnan_ctLR(i) ~= 0
%         c=CTLR_2(firstnan_ctLR(i):lastnan_ctLR(i),i);
%         
%         a=repmat(c(1),num,1);
%         b=repmat(c(end),num,1);
%         d=[a;c;b];
%         f=gaussfilt(1:length(d),d,win);
%         CTLR_3(firstnan_ctLR(i):lastnan_ctLR(i),i)=f(num+1:end-num);
%     end
end

DensHR_3=gsw_sigma0(SAHR_3,CTHR_3);
% DensLR_3=gsw_sigma0(SALR_3,CTLR_3);
%% save in netcdf 
% copy du fichier netcdf et on remplace la salinite et la temperature par
% les nouvelles matrices
% pour les donnees hr la matrice est de meme taille
if hr==1
    copyfile([conf.datadir_prof_hr list_tag(ll).name(1:end-8) '_CTDHR_prof.nc'],[conf.datadir_prof_sal_cor_hr list_tag(ll).name(1:end-8) '_sal_cor_CTDHR_prof.nc']);
    ficout=[conf.datadir_prof_sal_cor_hr list_tag(ll).name(1:end-8) '_sal_cor_CTDHR_prof.nc'];
    ncwrite(ficout,'TEMP_ADJUSTED', single(CTHR_3));
    ncwrite(ficout,'PSAL_ADJUSTED', single(SAHR_3));
    
else
    isfluo=0;
    isoxy=0;
    ficout=[conf.datadir,info_deployment.PI,'/',info_deployment.EXP '_effet_thermal_mass/' list_tag(ll).name(1:end-8) '_sal_cor_prof.nc'];
    ARGO_create(ficout,size(SAHR_2,2),length(Pi),isfluo,isoxy,0);
    %ncwrite(ficoutind,'PLATFORM_NUMBER',platform_number(:,I));
    ncwrite(ficout,'PLATFORM_NUMBER',MqHR.PLATFORM_NUMBER);
    ncwrite(ficout,'PI_NAME',MqHR.PI_NAME);
    ncwrite(ficout,'PROJECT_NAME',ncread(name_prof,'PROJECT_NAME'));
    ncwrite(ficout,'CYCLE_NUMBER',MqHR.CYCLE_NUMBER);
    Tqc=repmat('9',length(Pi),size(SAHR_2,2));
    Tqc(~isnan(Temp))='1';

    ncwrite(ficout,'JULD',MqHR.JULD);
    ncwrite(ficout,'JULD_LOCATION',MqHR.JULD_LOCATION);
    ncwrite(ficout,'LATITUDE', MqHR.LATITUDE);
    ncwrite(ficout,'LONGITUDE', MqHR.LONGITUDE);
    ncwrite(ficout,'PRES', single(Pres));
    ncwrite(ficout,'PRES_QC',Tqc);
    ncwrite(ficout,'TEMP', Temp);
    ncwrite(ficout,'TEMP_QC', Tqc);
    ncwrite(ficout,'PSAL', Sal);
    ncwrite(ficout,'PSAL_QC', Tqc);
    if isfluo
        ncwrite(ficout,'CHLA', Fluo(1:Nlevels,I));
        ncwrite(ficout,'CHLA_QC', Fqc(1:Nlevels,I));
    end
    if isoxy
        ncwrite(ficout,'DOXY', single(Oxy(1:Nlevels,I)));
        ncwrite(ficout,'DOXY_QC', Oqc(1:Nlevels,I));
    end
    ncwrite(ficout,'TEMP_ADJUSTED', single(CTHR_3));
    ncwrite(ficout,'PSAL_ADJUSTED', single(SAHR_3));
    ncwrite(ficout,'PRES_ADJUSTED', single(Pres));
    ncwrite(ficout,'PRES_ADJUSTED_QC',Tqc);
    ncwrite(ficout,'TEMP_ADJUSTED_QC', Tqc);
    ncwrite(ficout,'PSAL_ADJUSTED_QC', Tqc);

    ncwriteatt(ficout,'/','comment','TEMP : sc_correct_offset, TEMP_ADJUSTED : lia correction ');
    name_prof=ficout;
    smru_name = list_tag(ll).name(1:end-8);
    sc_write_global_attribute
end

% pour les donnees lr il faut re interpoller sur la grille lr avant
% d'enregistrer

% save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinit�/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_SALR_gf.mat'],'SALR_3')
% save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinit�/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_SAHR_gf.mat'],'SAHR_3')
% save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinit�/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_CTLR_gf.mat'],'CTLR_3')
% save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinit�/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_CTHR_gf.mat'],'CTHR_3')
% save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinit�/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_DensLR_gf.mat'],'DensLR_3')
% save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinit�/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_DensHR_gf.mat'],'DensHR_3')
