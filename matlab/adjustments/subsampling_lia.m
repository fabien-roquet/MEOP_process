%%%%%%%%%%%%%%Correct salinity and Temperature using sal_cor.m and Barker - July 2017

addpath(genpath('/Users/lsiegelma/Downloads/gsw_matlab_v3_05_7'));
addpath(genpath('/Users/lsiegelma/Documents/PhD/Matlab'));

clc
clf
clear all;
%close all;
%% PARAMETERS TO CHOOSE 

%plot tous les profils 1 ou 0
all_plot=1;
% time threshold to match profiles 0.0069 is 10 min %0.014is 20 min, 0.01 is 15 min, 0.021 is 30 min , 0.03125 % corresponds to 45min in julian day 
thres_min = 0.0069;
%%
format loose
format compact
dirplot='/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/plots/';

indivi = ['33';'35';'48';'49';'50'];
ind = 3;
individu = indivi(ind,:);

nc_nameHR = ['/Users/lsiegelma/Documents/PhD/Matlab/Data/TS_HR_netcdfs/ct112-0' individu '-14_CTDHR_prof.nc'];
nc_nameLR = ['/Users/lsiegelma/Documents/PhD/Matlab/Data/MEOP-CTD_2016-07-12/FRANCE/ncARGO/ct112-0' individu '-14_prof.nc'];
MqHR=ARGO_load_qc(nc_nameHR,1);
MqLR=ARGO_load_qc(nc_nameLR,1);
nc_att=ncloadatt_struct(nc_nameHR);

%Load High Resolution data
lat  = ncread(nc_nameHR,'LATITUDE');
lon  = ncread(nc_nameHR,'LONGITUDE');
Temp  = ncread(nc_nameHR,'TEMP_ADJUSTED');
Sal  = ncread(nc_nameHR,'PSAL_ADJUSTED');
Pres  = ncread(nc_nameHR,'PRES_ADJUSTED');
Date  = ncread(nc_nameHR,'JULD_LOCATION');

%Load Low Resolution data
latL  = ncread(nc_nameLR,'LATITUDE');
lonL  = ncread(nc_nameLR,'LONGITUDE');
TempL  = ncread(nc_nameLR,'TEMP_ADJUSTED');
SalL  = ncread(nc_nameLR,'PSAL_ADJUSTED');
PresL  = ncread(nc_nameLR,'PRES_ADJUSTED');
DateL  = ncread(nc_nameLR,'JULD_LOCATION');
%% quality check
%Make NaN missing values
latL(find(latL==99999))  = NaN;
lonL(find(lonL==99999))  = NaN;
TempL(find(TempL==99999))  = NaN;
SalL(find(SalL==99999))  = NaN;
PresL(find(PresL==99999))  = NaN;
DateL(find(DateL==99999))  = NaN;

% % remplace bad quality salinity values by NaN 
% psal_qc=ncread(nc_nameHR,'PSAL_ADJUSTED_QC');
% [x y]=find(psal_qc=='4');
% y=unique(y);

%% Match profiles with time threshold  (thres_min)
ind=1;
for k = 1:size(TempL,2)
    [m,DP]=min(abs((Date-DateL(k))));
    b(k)=m;
    if m<thres_min
        d_date(ind)=Date(DP);
        ind_hr(ind)=DP;
        ind_lr(ind)=k;
        ind=ind+1;
    end
end
%% check if some LR profiles are deeper than HR profiles 

A=PresL(:,ind_lr);

B = ~isnan(A);
% indices
Indices = arrayfun(@(x) find(B(:, x), 1, 'last'), 1:size(A, 2));
% values
lastnan_LR = arrayfun(@(x,y) A(x,y), Indices, 1:size(A, 2));

lastnan_tempHR=find_ndim(~isnan(Temp(:,ind_hr)),1,'last');
lastnan_salHR=find_ndim(~isnan(Sal(:,ind_hr)),1,'last');


ind_ct =  find((lastnan_LR>lastnan_tempHR));% & (lastnan_LR-lastnan_tempHR>10));
ind_sa =  find((lastnan_LR>lastnan_salHR));% & (lastnan_LR-lastnan_salHR>10));
ind_all = unique([ind_ct,ind_sa]); 

lon_all = lon(ind_hr(ind_all));
lonL_all = lonL(ind_lr(ind_all));

lat_all = lat(ind_hr(ind_all));
latL_all = latL(ind_lr(ind_all));

count = 1;
if size(ind_all,1)==0
    disp(['STEP ',num2str(count) ,': HR depth max > LR depth max for all profiles matched'])
end
while size(ind_all,1)~=0
    if (abs(lat_all-latL_all)>1e-4)  | (abs(lon_all-lonL_all)>1e-4)
        disp(['STEP ',num2str(count) , ' : for ',num2str(size(ind_all,2)), ' profiles LR depth max > HR depth max and some lon or lat differ.'])% Profiles # are : ', num2str(ind_all)])
    else
        disp(['STEP ' ,num2str(count) , ' : for ',num2str(size(ind_all,2)), ' profiles LR depth max > HR depth max but lon and lat always match exactly.'])% Profiles # are : ', num2str(ind_all)])
    end
    ind_hr(ind_all)=ind_hr(ind_all)-1;
    lastnan_tempHR=find_ndim(~isnan(Temp(:,ind_hr)),1,'last');
    lastnan_salHR=find_ndim(~isnan(Sal(:,ind_hr)),1,'last');
    ind_ct=[]; ind_sa=[]; ind_all=[];
    ind_ct =  find((lastnan_LR>lastnan_tempHR));% & (lastnan_LR-lastnan_tempHR>10));
    ind_sa =  find((lastnan_LR>lastnan_salHR));% & (lastnan_LR-lastnan_salHR>10));
    ind_all = unique([ind_ct;ind_sa]);
    count=count+1;
end
disp(['STEP ',num2str(count) ,': HR depth max > LR depth max for all profiles matched'])

%% removed wrong profile-pairs based on raw data quality check 

% remplace bad quality salinity values by NaN 
psal_qc=ncread(nc_nameHR,'PSAL_ADJUSTED_QC');
[x y]=find(psal_qc=='4');
y=unique(y);
y=intersect(y,ind_hr);
ind_rem=[];
for i=1:size(y,1)
    ind_rem(i)=find(ind_hr==y(i));
end

disp(['Numbers of pairs (LR+HR) of profiles removed : ',num2str(size(y,1)) , '/', num2str(size(ind_lr,2)), ', based on the HR PSAL quality flag (flag equals to 4)'])
ind_hr(ind_rem)=[];
ind_lr(ind_rem)=[];
d_date(ind_rem)=[];
%% plot non_matched profiles and their neighbors (has to be done for count = 1, i.e. before the while loop in the above section)
% for k=1:size(ind_all,2)
% clf;
% 
% subplot(1,2,1)
% plot(Temp(:,ind_hr(ind_all(k))-1:ind_hr(ind_all(k))-1),Pres(:,ind_hr(ind_all(k))),'k-')%,Temp(:,ind_hr(ind_all(k))+1:ind_hr(ind_all(k))+1),Pres(:,ind_hr(ind_all(k))),'k-')
% hold on 
% plot(Temp(:,ind_hr(ind_all(k))),Pres(:,ind_hr(ind_all(k))),'r-', 'linewidth',2)
% plot(TempL(:,ind_lr(ind_all(k))),PresL(:,ind_lr(ind_all(k))),'b-', 'linewidth',2)
% xlabel('HR Temperature (CT)')
% xlim([min(Temp(:,ind_hr(ind_all(k))))-1 max(Temp(:,ind_hr(ind_all(k))))+1])
% set(gca, 'YDir', 'reverse')
% grid on 
% 
% 
% subplot(1,2,2)
% plot(Sal(:,ind_hr(ind_all(k))-1:ind_hr(ind_all(k))-1),Pres(:,ind_hr(ind_all(k))),'k-')%,Sal(:,ind_hr(ind_all(k))+1:ind_hr(ind_all(k))+1),Pres(:,ind_hr(ind_all(k))),'k-')
% hold on 
% plot(Sal(:,ind_hr(ind_all(k))),Pres(:,ind_hr(ind_all(k))),'r-', 'linewidth',2)
% plot(SalL(:,ind_lr(ind_all(k))),PresL(:,ind_lr(ind_all(k))),'b-', 'linewidth',2)
% xlabel('HR Salinity (SA)')
% xlim([min(Sal(:,ind_hr(ind_all(k))))-0.5 max(Sal(:,ind_hr(ind_all(k))))+0.5])
% set(gca, 'YDir', 'reverse')
% grid on 
% 
% suptitle(['Profil # ',num2str(k),'/',num2str(size(ind_all,1)),'LR depth: ', num2str(lastnan_LR(ind_all(k))), ' Temp HR depth: ', num2str(lastnan_tempHR(ind_all(k))), ' Sal HR depth: ', num2str(lastnan_salHR(ind_all(k)))])
% pause(1)
% end

% clf
% plot(SalHR_brut(:,405),Pres(:,ind_hr(405)),'k-')%,Sal(:,ind_hr(ind_all(k))+1:ind_hr(ind_all(k))+1),Pres(:,ind_hr(ind_all(k))),'k-')
% hold on 
% plot(SalLR_brut(:,405),Pres(:,ind_hr(405)),'b-') 
% 
% xlabel('HR Salinity (SA)')
% %xlim([min(Sal(:,ind_hr(ind_all(k))))-0.5 max(Sal(:,ind_hr(ind_all(k))))+0.5])
% set(gca, 'YDir', 'reverse')
% grid on

%% apply step 1 to selected pairs of matching profiles: Thermal Mass Effect - Sal_cor with alpha=beta=0.03 (Vigan's coefficients)

TempHR_cor=[];
SalHR_cor=[];
DensHR_cor=[];
TempHR_brut=[];
SalHR_brut=[];
DensHR_brut=[];
SAHR_cor=[];
CTHR_cor=[];
SAHR_brut=[];
CTHR_brut=[];

TempLR_brut=[];
TempLR_cor=[];
SalLR_cor=[];
SalLR_brut=[];
DensLR_cor=[];
DensLR_brut=[];
SALR_cor=[];
CTLR_cor=[];
SALR_brut=[];
CTLR_brut=[];

nbpt_brut=[];
d_date=[];

alpha=0.03;
beta=0.03;
SI=1;
a=(2*alpha)./((SI*beta)+2);
b=1-((2*a)/alpha);
disp(['TME parameters: SI=1, alpha = beta = 0.03, a = ',num2str(a),' and b = ', num2str(b)])

for k = 1:size(ind_hr,2)
        %% correction HR
        T = Temp(:,ind_hr(k));
        S = Sal(:,ind_hr(k));
        p = Pres(:,ind_hr(k));
        d_date(k)=Date(ind_hr(k));
         
        saHR_b=gsw_SA_from_SP(S,p,lon(ind_hr(k)),lat(ind_hr(k))); % Absolute salinity (SA) from practical salinity (SP)
        ctHR_b=gsw_CT_from_t(saHR_b,T,p); % Conservative Temperature from in-situ temperature
        D=gsw_sigma0(saHR_b,ctHR_b);
        
        [Sc Tc] = sal_cor_20161214(alpha,beta,SI,T,S,p);
 
        sa=gsw_SA_from_SP(Sc,p,lon(ind_hr(k)),lat(ind_hr(k))); % Absolute salinity (SA) from practical salinity (SP)
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
        
        %% correction LR
        TL = TempL(:,ind_lr(k));
        SL = SalL(:,ind_lr(k));
        pL = PresL(:,ind_lr(k));
       % ind_lr(ind)=k;
        
%         sa_trans = gsw_SA_from_SP(SL,pL,lonL(k),latL(k));
%         ct_trans = gsw_CT_from_t(sa_trans,TL,pL);
%         dens_trans = gsw_sigma0(sa_trans,ct_trans);

        I = find(~isnan(pL));
        nbpt_brut(k,1)=length(I);
        Pi=(0:1000)';
        TiL=interp1(pL(I),TL(I),Pi);
        SiL=interp1(pL(I),SL(I),Pi);

        
        saLR_b=gsw_SA_from_SP(SiL,Pi,lonL(ind_lr(k)),latL(ind_lr(k))); % Absolute salinity (SA) from practical salinity (SP)
        %saLR_b = interp1(pL(I),sa_trans(I),Pi);
        ctLR_b=gsw_CT_from_t(saLR_b,TiL,Pi); % Conservative Temperature from in-situ temperature
        DL=gsw_sigma0(saLR_b,ctLR_b);
        %DL = interp1(pL(I),dens_trans(I),Pi);
        
        [ScL TcL] = sal_cor_20161214(alpha,beta,SI,TiL,SiL,Pi);

        saL=gsw_SA_from_SP(ScL,Pi,lonL(ind_lr(k)),latL(ind_lr(k))); % Absolute salinity (SA) from practical salinity (SP)
        ctL=gsw_CT_from_t(saL,TiL,Pi); % Conservative Temperature from in-situ temperature
        DcL=gsw_sigma0(saL,ctL);

        SALR_cor(:,k)=saL;
        CTLR_cor(:,k)=ctL;
        TempLR_cor(:,k)=TcL;
        SalLR_cor(:,k)=ScL;
        DensLR_cor(:,k)=DcL;
        SALR_brut(:,k)=saLR_b;
        CTLR_brut(:,k)=ctLR_b;
        TempLR_brut(:,k)=TiL;
        SalLR_brut(:,k)=SiL;
        DensLR_brut(:,k)=DL;
        %% Plotting
        % plot profil HR et LR avec TME
     
        if (all_plot & mod(k,5)==0   & sum(~isnan(ctHR_b))>0)
        mkdir([dirplot nc_att.smru_platform_code '/profils/step1/']);
        subplot(2,3,1)
        plot(ctHR_b,p,'k-',ct,p,'r-')
        xlabel('HR Temperature (CT)')
        xlim([min(ctHR_b)-0.5 max(ctHR_b)+0.5])
        grid on 
        set(gca, 'YDir', 'reverse')
        ylim([0 1000])
        subplot(2,3,2)
        plot(saHR_b,p,'k-',sa,p,'r-')
        xlabel('HR Salinity (SA)')
        xlim([min(saHR_b)-0.3 max(saHR_b)+0.3])
        set(gca, 'YDir', 'reverse')
        ylim([0 1000])
        grid on 
        subplot(2,3,3)
        plot(D,p,'k-',Dc,p,'r-')
        xlabel('HR Sigma 0')
        xlim([min(D)-0.3 max(D)+0.3])
        set(gca, 'YDir', 'reverse')
        ylim([0 1000])
        grid on 
        legend('raw','step 1','Location','southwest')

        subplot(2,3,4)
        plot(ctLR_b,p,'k-',ctL,p,'r-')
        xlabel('LR Temperature (CT)')
        xlim([min(ctHR_b)-0.5 max(ctHR_b)+0.5])
        set(gca, 'YDir', 'reverse')
        ylim([0 1000])
        grid on 
        subplot(2,3,5)
        plot(saLR_b,p,'k-',saL,p,'r-')
        xlabel('LR Salinity (SA)')
        xlim([min(saHR_b)-0.3 max(saHR_b)+0.3])
        set(gca, 'YDir', 'reverse')
        ylim([0 1000])
        grid on 
        subplot(2,3,6)
        plot(DL,p,'k-',DcL,p,'r-')
        xlabel('LR Sigma 0')
        xlim([min(D)-0.3 max(D)+0.3])
        set(gca, 'YDir', 'reverse')
        grid on 
        ylim([0 1000])

        print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/profils/step1/' nc_att.smru_platform_code '_profilHR_' num2str(k)]);
        end 
end
clear saHR_b ctHR_b D T S Sc p sa ct Dc saLR_b ctLR_b TL SL ScL SiL pL saL ctL DcL DP k m ind I Tc TcL TiL DL A B count

%% Step #2 - Density Inversion Removal (Paul Barker and Trevor McDougall, 2017)
SAHR_2=gsw_stabilise_SA_const_t(SAHR_cor,CTHR_cor,Pres(:,ind_hr),lon(ind_hr),lat(ind_hr));
CTHR_2=gsw_CT_from_t(SAHR_2,TempHR_brut,Pres(:,ind_hr));
DensHR_2=gsw_sigma0(SAHR_2,CTHR_2);

SALR_2=gsw_stabilise_SA_const_t(SALR_cor,CTLR_cor,Pres(:,ind_lr),lonL(ind_lr),latL(ind_lr));
CTLR_2=gsw_CT_from_t(SALR_2,TempLR_brut,Pres(:,ind_lr));
DensLR_2=gsw_sigma0(SALR_2,CTLR_2);
%% checking plots
if all_plot
mkdir([dirplot nc_att.smru_platform_code '/profils/step2/']);
for ind=1:length(ind_hr)
    if (all_plot & mod(ind,5)==0)

        % find prof_max associated to selected profil
        validIndices = ~isnan(TempHR_brut(:,ind));
        v=+validIndices;  
        y = v; y(v == 0) = [];
        last = y(end);
        [r,c] =find(v==last);
        pmax=max(r);

        % plot T, S, Rho corrected and raw

        p=Pres(:,ind);
        subplot(2,3,1)
            plot(CTHR_brut(:,ind),p,'k-',CTHR_cor(:,ind),p,'r-',CTHR_2(:,ind),p,'b-')
            xlabel('CT (°C)')
            set(gca, 'YDir', 'reverse')
            grid on
            ylim([0 pmax])
            x1=get(gca,'xlim');
        subplot(2,3,2)
            plot(SAHR_brut(:,ind),p,'k-',SAHR_cor(:,ind),p,'r-',SAHR_2(:,ind),p,'b-')
            xlabel('SA (g.kg^{-1})')
            set(gca, 'YDir', 'reverse')
            grid on
            ylim([0 pmax])
            x2=get(gca,'xlim');
            title(['HR Profil #',num2str(ind)],'FontSize',18)
        subplot(2,3,3)       
            plot(DensHR_brut(:,ind),p,'k-',DensHR_cor(:,ind),p,'r-',DensHR_2(:,ind),p,'b-')
            xlabel('\rho_{0} (kg.m^{-3})')
            set(gca, 'YDir', 'reverse')
            grid on
            ylim([0 pmax])
            x3=get(gca,'xlim');
            legend({'raw data','sal cor','sal cor + Mc Dougal'},'Location','southwest','fontsize',5)
        subplot(2,3,4)
            plot(CTLR_brut(:,ind),p,'k-',CTLR_cor(:,ind),p,'r-',CTLR_2(:,ind),p,'b-')
            xlabel('CT (°C)')
            set(gca, 'YDir', 'reverse')
            grid on
            ylim([0 pmax])
            xlim(x1)
        subplot(2,3,5)
            plot(SALR_brut(:,ind),p,'k-',SALR_cor(:,ind),p,'r-',SALR_2(:,ind),p,'b-')
            xlabel('SA (g.kg^{-1})')
            set(gca, 'YDir', 'reverse')
            grid on
            ylim([0 pmax])
            xlim(x2)
            title(['LR'],'FontSize',18)
        subplot(2,3,6)       
            plot(DensLR_brut(:,ind),p,'k-',DensLR_cor(:,ind),p,'r-',DensLR_2(:,ind),p,'b-')
            xlabel('\rho_{0} (kg.m^{-3})')
            set(gca, 'YDir', 'reverse')
            grid on
            ylim([0 pmax])
            xlim(x3)
        print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/profils/step2/' nc_att.smru_platform_code '_profilHR_' num2str(ind)]);   
    end
end
end
%% save variables raw, step 1, step 2
mkdir(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code])
% save Temp & Sal brut, HR & LR
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_TempLR_raw.mat'],'TempLR_brut')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_TempHR_raw.mat'],'TempHR_brut')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_SalLR_raw.mat'],'SalLR_brut')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_SalHR_raw.mat'],'SalHR_brut')
% save CT,SA & Sigma_0 brut, HR & LR
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_SALR_raw.mat'],'SALR_brut')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_SAHR_raw.mat'],'SAHR_brut')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_CTLR_raw.mat'],'CTLR_brut')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_CTHR_raw.mat'],'CTHR_brut')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_DensLR_raw.mat'],'DensLR_brut')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_DensHR_raw.mat'],'DensHR_brut')
% save CT,SA & Sigma_0 TME, HR & LR
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_SALR_tme.mat'],'SALR_cor')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_SAHR_tme.mat'],'SAHR_cor')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_CTLR_tme.mat'],'CTLR_cor')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_CTHR_tme.mat'],'CTHR_cor')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_DensLR_tme.mat'],'DensLR_cor')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_DensHR_tme.mat'],'DensHR_cor')
% save CT,SA & Sigma_0 DIR, HR & LR
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_SALR_dir.mat'],'SALR_2')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_SAHR_dir.mat'],'SAHR_2')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_CTLR_dir.mat'],'CTLR_2')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_CTHR_dir.mat'],'CTHR_2')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_DensLR_dir.mat'],'DensLR_2')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_DensHR_dir.mat'],'DensHR_2')
%% Step #3 - Gaussian filter smoothing 
firstnan_ctHR=find_ndim(~isnan(CTHR_2),1,'first');
lastnan_ctHR=find_ndim(~isnan(CTHR_2),1,'last');

firstnan_saHR=find_ndim(~isnan(SAHR_2),1,'first');
lastnan_saHR=find_ndim(~isnan(SAHR_2),1,'last');


firstnan_ctLR=find_ndim(~isnan(CTLR_2),1,'first');
lastnan_ctLR=find_ndim(~isnan(CTLR_2),1,'last');

firstnan_saLR=find_ndim(~isnan(SALR_2),1,'first');
lastnan_saLR=find_ndim(~isnan(SALR_2),1,'last');

lastnan_densHR=find_ndim(~isnan(DensHR_2),1,'last');
%% find if a nan value is in between non nan-values and set to nan all values in between NaNs so that the gaussian filter doesnt propagate NaNs more than necessary
for i=1:size(SAHR_2,2)
    if sum(isnan(SAHR_2(firstnan_saHR(i):lastnan_saHR(i),i)))~= 0
        x=find(isnan(SAHR_2(firstnan_saHR(i):lastnan_saHR(i),i)));
        row=max(x);
        SAHR_2(1:x,i)=NaN;
    end
end
for i=1:size(CTHR_2,2)
    if sum(isnan(CTHR_2(firstnan_ctHR(i):lastnan_ctHR(i),i)))~= 0
        x=find(isnan(CTHR_2(firstnan_ctHR(i):lastnan_ctHR(i),i)));
        row=max(x);
        CTHR_2(1:x,i)=NaN;
    end
end
for i=1:size(DensHR_2,2)
    if sum(isnan(SAHR_2(firstnan_saHR(i):lastnan_saHR(i),i)))~= 0
        x=find(isnan(SAHR_2(firstnan_saHR(i):lastnan_saHR(i),i)));
        row=max(x);
        DensHR_2(1:x,i)=NaN;
    end
end
firstnan_saHR=find_ndim(~isnan(SAHR_2),1,'first');
firstnan_ctHR=find_ndim(~isnan(CTHR_2),1,'first');
firstnan_densHR=find_ndim(~isnan(DensHR_2),1,'first');
%% ordre 1 
SAHR_3=SAHR_2*NaN;
CTHR_3=CTHR_2*NaN;

SALR_3=SALR_2*NaN;
CTLR_3=CTLR_2*NaN;

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
    
    if firstnan_saLR(i) ~= 0
        c=SALR_2(firstnan_saLR(i):lastnan_saLR(i),i);
        
        a=repmat(c(1),num,1);
        b=repmat(c(end),num,1);
        d=[a;c;b];
        s=gaussfilt(1:length(d),d,win);
        SALR_3(firstnan_saLR(i):lastnan_saLR(i),i)=s(num+1:end-num);
    end
        
    if firstnan_ctLR(i) ~= 0
        c=CTLR_2(firstnan_ctLR(i):lastnan_ctLR(i),i);
        
        a=repmat(c(1),num,1);
        b=repmat(c(end),num,1);
        d=[a;c;b];
        f=gaussfilt(1:length(d),d,win);
        CTLR_3(firstnan_ctLR(i):lastnan_ctLR(i),i)=f(num+1:end-num);
    end
end

DensHR_3=gsw_sigma0(SAHR_3,CTHR_3);
DensLR_3=gsw_sigma0(SALR_3,CTLR_3);
%% save variables step 3
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_SALR_gf.mat'],'SALR_3')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_SAHR_gf.mat'],'SAHR_3')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_CTLR_gf.mat'],'CTLR_3')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_CTHR_gf.mat'],'CTHR_3')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_DensLR_gf.mat'],'DensLR_3')
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_DensHR_gf.mat'],'DensHR_3')

%% checking plot

if all_plot
mkdir([dirplot nc_att.smru_platform_code '/profils/step3/']);
for k=1:length(ind_hr)
   p=Pres(:,k);
   pmax=lastnan_densHR(k);pmin=0;
   if (all_plot & mod(k,5)==0 & pmax>60)
       clf;
    % plot T, S, Rho corrected and raw
    subplot(2,3,1)
        plot(CTHR_brut(:,k),p,'k-',CTHR_cor(:,k),p,'g',CTHR_2(:,k),p,'b-',CTHR_3(:,k),p,'r')
        xlabel('CT (°C)')
        set(gca, 'YDir', 'reverse')
        grid on
        ylim([pmin pmax])
        x1=get(gca,'xlim');
    subplot(2,3,2)
        plot(SAHR_brut(:,k),p,'k-',SAHR_cor(:,k),p,'g',SAHR_2(:,k),p,'b-',SAHR_3(:,k),p,'r-')
        xlabel('SA (g.kg^{-1})')
        set(gca, 'YDir', 'reverse')
        grid on
        ylim([pmin pmax])
        legend({'raw','step 1','step 2','step 3'},'Location','southwest', 'fontsize',5)
        title(['HR Profil #',num2str(k) ],'FontSize',14)
        x2=get(gca,'xlim');
    subplot(2,3,3)       
        plot(DensHR_brut(:,k),p,'k-',DensHR_cor(:,k),p,'g',DensHR_2(:,k),p,'b-',DensHR_3(:,k),p,'r-')
        xlabel('\rho_{0} (kg.m^{-3})')
        set(gca, 'YDir', 'reverse')
        grid on
        ylim([pmin pmax])
        x3=get(gca,'xlim');
    subplot(2,3,4)
        plot(CTLR_brut(:,k),p,'k-',CTLR_cor(:,k),p,'g',CTLR_2(:,k),p,'b-',CTLR_3(:,k),p,'r')
        xlabel('CT (°C)')
        set(gca, 'YDir', 'reverse')
        grid on
        ylim([pmin pmax])
        xlim(x1)
    subplot(2,3,5)
        plot(SALR_brut(:,k),p,'k-',SALR_cor(:,k),p,'g',SALR_2(:,k),p,'b-',SALR_3(:,k),p,'r-')
        xlabel('SA (g.kg^{-1})')
        set(gca, 'YDir', 'reverse')
        grid on
        ylim([pmin pmax])
        xlim(x2)
        title('LR','FontSize',14)
    subplot(2,3,6)       
        plot(DensLR_brut(:,k),p,'k-',DensLR_cor(:,k),p,'g',DensLR_2(:,k),p,'b-',DensLR_3(:,k),p,'r-')
        xlabel('\rho_{0} (kg.m^{-3})')
        set(gca, 'YDir', 'reverse')
        grid on
        ylim([pmin pmax])
        xlim(x3)
        print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/profils/step3/' nc_att.smru_platform_code '_profilHR_' num2str(k)]);
   end
end
end
%% Comparison HR LR RAW DATA

if all_plot
    mkdir([dirplot nc_att.smru_platform_code '/comparison/raw/']);
    for i=1:10:size(CTHR_brut,2)
        %% figure HR raw/cor et smooth
        clf
        subplot(2,3,1)
        plot(CTHR_brut(:,i),p,'g-',CTLR_brut(:,i),p,'b-')
        xlabel('CT (°C)')
        ylabel('Depth (m)')
        if sum(~isnan(CTHR_brut(:,i)))
            xlim([min(CTHR_brut(:,i))-0.5 max(CTHR_brut(:,i))+0.5])
        else
            xlim([-0.5 5])
        end
        
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,2)
        plot(SAHR_brut(:,i),p,'g-',SALR_brut(:,i),p,'b-')
        xlabel('SA (g.kg^{-1})')
        if sum(~isnan(SAHR_brut(:,i)))
            xlim([min(SAHR_brut(:,i))-0.1 max(SAHR_brut(:,i))+0.1])
        else
            xlim([33 35])
        end
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,3)
        plot(DensHR_brut(:,i),p,'g-',DensLR_brut(:,i),p,'b-')
        xlabel('\rho_{0} (kg.m^{-3})')
        if sum(~isnan(DensHR_brut(:,i)))
            xlim([min(DensHR_brut(:,i))-0.1 max(DensHR_brut(:,i))+0.1])
        else
        xlim([26 27.5])
        end
        legend('HR','LR','Location','southwest')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,4)
        plot(CTHR_brut(:,i)-CTLR_brut(:,i),p,'k-')
        T=find(~isnan(CTHR_brut(:,i)) & ~isnan(CTLR_brut(:,i)));
        ecart_T(i,1)=rms(CTHR_brut(T,i)-CTLR_brut(T,i));
        anno=['rms : ',num2str(ecart_T(i,1))];
        title(anno);
        xlabel('CT HR - CT LR')
        ylabel('Depth (m)')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,5)
        plot(SAHR_brut(:,i)-SALR_brut(:,i),p,'k-')
        T=find(~isnan(SAHR_brut(:,i)) & ~isnan(SALR_brut(:,i)));
        ecart_S(i,1)=rms(SAHR_brut(T,i)-SALR_brut(T,i));
        anno=['rms : ',num2str(ecart_S(i,1))];
        title(anno);
        xlabel('SA HR - SA LR')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,6)
        plot(DensHR_brut(:,i)-DensLR_brut(:,i),p,'k-')
        T=find(~isnan(DensHR_brut(:,i)) & ~isnan(DensLR_brut(:,i)));
        ecart_D(i,1)=rms(DensHR_brut(T,i)-DensLR_brut(T,i));
        anno=['rms : ',num2str(ecart_D(i,1))];
        title(anno);
        xlabel('Dens HR - Dens LR')
        set(gca, 'YDir', 'reverse')
        grid on
        
        suptitle(['RAW DATA - ' nc_att.smru_platform_code ' - Profil # '  num2str(i) ' - # points LR: ' num2str(nbpt_brut(i))]);
        print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/raw/' nc_att.smru_platform_code '_profil_' num2str(i)]);
    end
end
%% Comparison HR LR STEP 1 THERMAL MASS EFFECT
if all_plot
    mkdir([dirplot nc_att.smru_platform_code '/comparison/step1/']);
    for i=1:10:size(CTHR_cor,2)
        %% figure HR raw/cor et smooth
        clf
        subplot(2,3,1)
        plot(CTHR_cor(:,i),p,'g-',CTLR_cor(:,i),p,'b-')
        xlabel('CT (°C)')
        ylabel('Depth (m)')
        if sum(~isnan(CTHR_cor(:,i)))
            xlim([min(CTHR_cor(:,i))-0.5 max(CTHR_cor(:,i))+0.5])
        else
            xlim([-0.5 5])
        end
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,2)
        plot(SAHR_cor(:,i),p,'g-',SALR_cor(:,i),p,'b-')
        xlabel('SA (g.kg^{-1})')
        if sum(~isnan(SAHR_cor(:,i)))
            xlim([min(SAHR_cor(:,i))-0.1 max(SAHR_cor(:,i))+0.1])
        else
            xlim([33 35])
        end
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,3)
        plot(DensHR_cor(:,i),p,'g-',DensLR_cor(:,i),p,'b-')
        xlabel('\rho_{0} (kg.m^{-3})')
        if sum(~isnan(DensHR_cor(:,i)))
            xlim([min(DensHR_cor(:,i))-0.1 max(DensHR_cor(:,i))+0.1])
        else
        xlim([26 27.5])
        end
        legend('HR TME','LR TME','Location','southwest')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,4)
        plot(CTHR_cor(:,i)-CTLR_cor(:,i),p,'k-')
        T=find(~isnan(CTHR_cor(:,i)) & ~isnan(CTLR_cor(:,i)));
        ecart_T(i,1)=rms(CTHR_cor(T,i)-CTLR_cor(T,i));
        anno=['rms : ',num2str(ecart_T(i,1))];
        title(anno);
        xlabel('CT HR - CT LR')
        ylabel('Depth (m)')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,5)
        plot(SAHR_cor(:,i)-SALR_cor(:,i),p,'k-')
        T=find(~isnan(SAHR_cor(:,i)) & ~isnan(SALR_cor(:,i)));
        ecart_S(i,1)=rms(SAHR_cor(T,i)-SALR_cor(T,i));
        anno=['rms : ',num2str(ecart_S(i,1))];
        title(anno);
        xlabel('SA HR - SA LR')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,6)
        plot(DensHR_cor(:,i)-DensLR_cor(:,i),p,'k-')
        T=find(~isnan(DensHR_cor(:,i)) & ~isnan(DensLR_cor(:,i)));
        ecart_D(i,1)=rms(DensHR_cor(T,i)-DensLR_cor(T,i));
        anno=['rms : ',num2str(ecart_D(i,1))];
        title(anno);
        xlabel('Dens HR - Dens LR')
        set(gca, 'YDir', 'reverse')
        grid on
        
        suptitle(['STEP 1 : THERMAL MASS EFFECT (TME) - ' nc_att.smru_platform_code ' - Profil # '  num2str(i) ' - # points LR: ' num2str(nbpt_brut(i))]);
        print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/step1/' nc_att.smru_platform_code '_profil_' num2str(i)]);
    end
end
%% Comparison HR LR STEP 2 DENSITY INVERSION REMOVAL

if all_plot
    mkdir([dirplot nc_att.smru_platform_code '/comparison/step2/']);
    for i=1:20:size(CTHR_2,2)
        %% figure HR raw/cor et smooth
        clf
        subplot(2,3,1)
        plot(CTHR_2(:,i),p,'g-',CTLR_2(:,i),p,'b-')
        xlabel('CT (°C)')
        ylabel('Depth (m)')
        if sum(~isnan(CTHR_2(:,i)))
            xlim([min(CTHR_2(:,i))-0.5 max(CTHR_2(:,i))+0.5])
        else
            xlim([-0.5 5])
        end
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,2)
        plot(SAHR_2(:,i),p,'g-',SALR_2(:,i),p,'b-')
        xlabel('SA (g.kg^{-1})')
        if sum(~isnan(SAHR_2(:,i)))
            xlim([min(SAHR_2(:,i))-0.1 max(SAHR_2(:,i))+0.1])
        else
            xlim([33 35])
        end
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,3)
        plot(DensHR_2(:,i),p,'g-',DensLR_2(:,i),p,'b-')
        xlabel('\rho_{0} (kg.m^{-3})')
        if sum(~isnan(DensHR_2(:,i)))
            xlim([min(DensHR_2(:,i))-0.1 max(DensHR_2(:,i))+0.1])
        else
        xlim([26 27.5])
        end
        legend('HR DIR','LR DIR','Location','southwest')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,4)
        plot(CTHR_2(:,i)-CTLR_2(:,i),p,'k-')
        T=find(~isnan(CTHR_2(:,i)) & ~isnan(CTLR_2(:,i)));
        ecart_T(i,1)=rms(CTHR_2(T,i)-CTLR_2(T,i));
        anno=['rms : ',num2str(ecart_T(i,1))];
        title(anno);
        xlabel('CT HR - CT LR')
        ylabel('Depth (m)')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,5)
        plot(SAHR_2(:,i)-SALR_2(:,i),p,'k-')
        T=find(~isnan(SAHR_2(:,i)) & ~isnan(SALR_2(:,i)));
        ecart_S(i,1)=rms(SAHR_2(T,i)-SALR_2(T,i));
        anno=['rms : ',num2str(ecart_S(i,1))];
        title(anno);
        xlabel('SA HR - SA LR')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,6)
        plot(DensHR_2(:,i)-DensLR_2(:,i),p,'k-')
        T=find(~isnan(DensHR_2(:,i)) & ~isnan(DensLR_2(:,i)));
        ecart_D(i,1)=rms(DensHR_2(T,i)-DensLR_2(T,i));
        anno=['rms : ',num2str(ecart_D(i,1))];
        title(anno);
        xlabel('Dens HR - Dens LR')
        set(gca, 'YDir', 'reverse')
        grid on
        
        suptitle(['STEP 2 : DENSITY INVERSION REMOVAL (DIR) - ' nc_att.smru_platform_code ' - Profil # '  num2str(i) ' - # points LR: ' num2str(nbpt_brut(i))]);
        print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/step2/' nc_att.smru_platform_code '_profil_' num2str(i)]);
    end
end
%% Comparison HR LR STEP 3 GAUSSIAN FILTER (GF)

if all_plot
    mkdir([dirplot nc_att.smru_platform_code '/comparison/step3/']);
    for i=1:20:size(CTHR_3,2)
        %% figure HR raw/cor et smooth
        clf
        subplot(2,3,1)
        plot(CTHR_3(:,i),p,'g-',CTLR_3(:,i),p,'b-')
        xlabel('CT (°C)')
        ylabel('Depth (m)')
        if sum(~isnan(CTHR_3(:,i)))
            xlim([min(CTHR_3(:,i))-0.5 max(CTHR_3(:,i))+0.5])
        else
            xlim([-0.5 5])
        end
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,2)
        plot(SAHR_3(:,i),p,'g-',SALR_3(:,i),p,'b-')
        xlabel('SA (g.kg^{-1})')
        if sum(~isnan(SAHR_3(:,i)))
            xlim([min(SAHR_3(:,i))-0.1 max(SAHR_3(:,i))+0.1])
        else
            xlim([33 35])
        end
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,3)
        plot(DensHR_3(:,i),p,'g-',DensLR_3(:,i),p,'b-')
        xlabel('\rho_{0} (kg.m^{-3})')
        if sum(~isnan(DensHR_3(:,i)))
            xlim([min(DensHR_3(:,i))-0.1 max(DensHR_3(:,i))+0.1])
        else
        xlim([26 27.5])
        end
        legend('HR GF','LR GF','Location','southwest')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,4)
        plot(CTHR_3(:,i)-CTLR_3(:,i),p,'k-')
        T=find(~isnan(CTHR_3(:,i)) & ~isnan(CTLR_3(:,i)));
        ecart_T(i,1)=rms(CTHR_3(T,i)-CTLR_3(T,i));
        anno=['rms : ',num2str(ecart_T(i,1))];
        title(anno);
        xlabel('CT HR - CT LR')
        ylabel('Depth (m)')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,5)
        plot(SAHR_3(:,i)-SALR_3(:,i),p,'k-')
        T=find(~isnan(SAHR_3(:,i)) & ~isnan(SALR_3(:,i)));
        ecart_S(i,1)=rms(SAHR_3(T,i)-SALR_3(T,i));
        anno=['rms : ',num2str(ecart_S(i,1))];
        title(anno);
        xlabel('SA HR - SA LR')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,6)
        plot(DensHR_3(:,i)-DensLR_3(:,i),p,'k-')
        T=find(~isnan(DensHR_3(:,i)) & ~isnan(DensLR_3(:,i)));
        ecart_D(i,1)=rms(DensHR_3(T,i)-DensLR_3(T,i));
        anno=['rms : ',num2str(ecart_D(i,1))];
        title(anno);
        xlabel('Dens HR - Dens LR')
        set(gca, 'YDir', 'reverse')
        grid on
        
        suptitle(['STEP 3 : GAUSSIAN FILTER (GF) - ' nc_att.smru_platform_code ' - Profil # '  num2str(i) ' - # points LR: ' num2str(nbpt_brut(i))]);
        print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/step3/' nc_att.smru_platform_code '_profil_' num2str(i)]);
    end
end
%% Comparison HR et LR RAW vs STEP1 (TME)

if all_plot
    mkdir([dirplot nc_att.smru_platform_code '/comparison/HR_LR/step1/']);
    for i=1:20:size(CTHR_cor,2)
        clf
        p=Pres(:,i);
        subplot(2,3,1)
        plot(CTHR_brut(:,i),p,'k-',CTHR_cor(:,i),p,'r-',CTLR_brut(:,i),p,'b-',CTLR_cor(:,i),p,'g-')
        xlabel('CT (°C)')
        ylabel('Depth (m)')
        if sum(~isnan(CTHR_cor(:,i)))
            xlim([min(CTHR_cor(:,i))-0.5 max(CTHR_cor(:,i))+0.5])
        else
            xlim([-0.5 5])
        end
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,2)
        plot(SAHR_brut(:,i),p,'k-',SAHR_cor(:,i),p,'r-',SALR_brut(:,i),p,'b-',SALR_cor(:,i),p,'g-')
        xlabel('SA (g.kg^{-1})')
        if sum(~isnan(SAHR_cor(:,i)))
            xlim([min(SAHR_cor(:,i))-0.1 max(SAHR_cor(:,i))+0.1])
        else
            xlim([33 35])
        end
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,3)
        plot(DensHR_brut(:,i),p,'k-',DensHR_cor(:,i),p,'r-',DensLR_brut(:,i),p,'b-',DensLR_cor(:,i),p,'g-')
        xlabel('\rho_{0} (kg.m^{-3})')
        if sum(~isnan(DensHR_cor(:,i)))
            xlim([min(DensHR_cor(:,i))-0.1 max(DensHR_cor(:,i))+0.1])
        else
        xlim([26 27.5])
        end
        legend('HR RAW','HR TME','LR RAW','LR TME','Location','southwest')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,4)
        plot(CTHR_brut(:,i)-CTHR_cor(:,i),p,'b-',CTLR_brut(:,i)-CTLR_cor(:,i),p,'r-')
        T=find(~isnan(CTHR_brut(:,i)) & ~isnan(CTHR_cor(:,i)));
        ecart_T(i,1)=rms(CTHR_brut(T,i)-CTHR_cor(T,i));
        T=find(~isnan(CTLR_brut(:,i)) & ~isnan(CTLR_cor(:,i)));
        ecart_L(i,1)=rms(CTLR_brut(T,i)-CTLR_cor(T,i));
        anno=['rms HR: ',num2str(ecart_T(i,1),'%.4f') ' rms LR: ' num2str(ecart_L(i,1),'%.4f') ];
        title(anno);
        xlabel('CT RAW - CT TME')
        ylabel('Depth (m)')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,5)
        plot(SAHR_brut(:,i)-SAHR_cor(:,i),p,'b-',SALR_brut(:,i)-SALR_cor(:,i),p,'r-')
        T=find(~isnan(SAHR_brut(:,i)) & ~isnan(SAHR_cor(:,i)));
        ecart_T(i,1)=rms(SAHR_brut(T,i)-SAHR_cor(T,i));
        T=find(~isnan(SALR_brut(:,i)) & ~isnan(SALR_cor(:,i)));
        ecart_L(i,1)=rms(SALR_brut(T,i)-SALR_cor(T,i));
        anno=['rms HR: ',num2str(ecart_T(i,1),'%.4f') ' rms LR: ' num2str(ecart_L(i,1),'%.4f') ];
        title(anno);
        xlabel('SA RAW - SA TME')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,6)
        plot(DensHR_brut(:,i)-DensHR_cor(:,i),p,'b-',DensLR_brut(:,i)-DensLR_cor(:,i),p,'r-')
        T=find(~isnan(DensHR_brut(:,i)) & ~isnan(DensHR_cor(:,i)));
        ecart_T(i,1)=rms(DensHR_brut(T,i)-DensHR_cor(T,i));
        T=find(~isnan(DensLR_brut(:,i)) & ~isnan(DensLR_cor(:,i)));
        ecart_L(i,1)=rms(DensLR_brut(T,i)-DensLR_cor(T,i));
        anno=['rms HR: ',num2str(ecart_T(i,1),'%.4f') ' rms LR: ' num2str(ecart_L(i,1),'%.4f') ];
        title(anno);
        xlabel('Dens RAW - Dens TME')
        set(gca, 'YDir', 'reverse')
        grid on
        legend('HR RAW - HR TME','LR RAW - LR TME','Location','southeast')
        
        suptitle(['HR & LR RAW VS STEP 1 (TME) - ' nc_att.smru_platform_code ' - Profile # '  num2str(i) ' - # points LR: ' num2str(nbpt_brut(i))]);
        set(gcf,'PaperPositionMode','auto')
        print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/HR_LR/step1/' nc_att.smru_platform_code '_profil_' num2str(i)]);
    end
end
%% Comparison HR et LR RAW vs STEP2 (DIR)

if all_plot
    mkdir([dirplot nc_att.smru_platform_code '/comparison/HR_LR/step2/']);
    for i=1:20:size(CTHR_2,2)
        clf
        p=Pres(:,i);
        subplot(2,3,1)
        plot(CTHR_brut(:,i),p,'k-',CTHR_2(:,i),p,'r-',CTLR_brut(:,i),p,'b-',CTLR_2(:,i),p,'g-')
        xlabel('CT (°C)')
        ylabel('Depth (m)')
        if sum(~isnan(CTHR_2(:,i)))
            xlim([min(CTHR_2(:,i))-0.5 max(CTHR_2(:,i))+0.5])
        else
            xlim([-0.5 5])
        end
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,2)
        plot(SAHR_brut(:,i),p,'k-',SAHR_2(:,i),p,'r-',SALR_brut(:,i),p,'b-',SALR_2(:,i),p,'g-')
        xlabel('SA (g.kg^{-1})')
        if sum(~isnan(SAHR_2(:,i)))
            xlim([min(SAHR_2(:,i))-0.1 max(SAHR_2(:,i))+0.1])
        else
            xlim([33 35])
        end
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,3)
        plot(DensHR_brut(:,i),p,'k-',DensHR_2(:,i),p,'r-',DensLR_brut(:,i),p,'b-',DensLR_2(:,i),p,'g-')
        xlabel('\rho_{0} (kg.m^{-3})')
        if sum(~isnan(DensHR_2(:,i)))
            xlim([min(DensHR_2(:,i))-0.1 max(DensHR_2(:,i))+0.1])
        else
        xlim([26 27.5])
        end
        legend('HR RAW','HR DIR','LR RAW','LR DIR','Location','southwest')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,4)
        plot(CTHR_brut(:,i)-CTHR_2(:,i),p,'b-',CTLR_brut(:,i)-CTLR_2(:,i),p,'r-')
        T=find(~isnan(CTHR_brut(:,i)) & ~isnan(CTHR_2(:,i)));
        ecart_T(i,1)=rms(CTHR_brut(T,i)-CTHR_2(T,i));
        T=find(~isnan(CTLR_brut(:,i)) & ~isnan(CTLR_2(:,i)));
        ecart_L(i,1)=rms(CTLR_brut(T,i)-CTLR_2(T,i));
        anno=['rms HR: ',num2str(ecart_T(i,1),'%.4f') ' rms LR: ' num2str(ecart_L(i,1),'%.4f') ];
        title(anno);
        xlabel('CT RAW - CT DIR')
        ylabel('Depth (m)')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,5)
        plot(SAHR_brut(:,i)-SAHR_2(:,i),p,'b-',SALR_brut(:,i)-SALR_2(:,i),p,'r-')
        T=find(~isnan(SAHR_brut(:,i)) & ~isnan(SAHR_2(:,i)));
        ecart_T(i,1)=rms(SAHR_brut(T,i)-SAHR_2(T,i));
        T=find(~isnan(SALR_brut(:,i)) & ~isnan(SALR_2(:,i)));
        ecart_L(i,1)=rms(SALR_brut(T,i)-SALR_2(T,i));
        anno=['rms HR: ',num2str(ecart_T(i,1),'%.4f') ' rms LR: ' num2str(ecart_L(i,1),'%.4f') ];
        title(anno);
        xlabel('SA RAW - SA DIR')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,6)
        plot(DensHR_brut(:,i)-DensHR_2(:,i),p,'b-',DensLR_brut(:,i)-DensLR_2(:,i),p,'r-')
        T=find(~isnan(DensHR_brut(:,i)) & ~isnan(DensHR_2(:,i)));
        ecart_T(i,1)=rms(DensHR_brut(T,i)-DensHR_2(T,i));
        T=find(~isnan(DensLR_brut(:,i)) & ~isnan(DensLR_2(:,i)));
        ecart_L(i,1)=rms(DensLR_brut(T,i)-DensLR_2(T,i));
        anno=['rms HR: ',num2str(ecart_T(i,1),'%.4f') ' rms LR: ' num2str(ecart_L(i,1),'%.4f') ];
        title(anno);
        xlabel('Dens RAW - Dens DIR')
        set(gca, 'YDir', 'reverse')
        grid on
        legend('HR RAW - HR DIR','LR RAW - LR DIR','Location','southeast')
        
        suptitle(['HR & LR RAW VS STEP 2 (DIR) - ' nc_att.smru_platform_code ' - Profile # '  num2str(i) ' - # points LR: ' num2str(nbpt_brut(i))]);
        set(gcf,'PaperPositionMode','auto')
        print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/HR_LR/step2/' nc_att.smru_platform_code '_profil_' num2str(i)]);
    end
end
%% Comparison HR et LR RAW vs STEP3 (GF)

if all_plot
    mkdir([dirplot nc_att.smru_platform_code '/comparison/HR_LR/step3/']);
    for i=1:20:size(CTHR_3,2)
        clf
        p=Pres(:,i);
        subplot(2,3,1)
        plot(CTHR_brut(:,i),p,'k-',CTHR_3(:,i),p,'r-',CTLR_brut(:,i),p,'b-',CTLR_3(:,i),p,'g-')
        xlabel('CT (°C)')
        ylabel('Depth (m)')
        if sum(~isnan(CTHR_3(:,i)))
            xlim([min(CTHR_3(:,i))-0.5 max(CTHR_3(:,i))+0.5])
        else
            xlim([-0.5 5])
        end
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,2)
        plot(SAHR_brut(:,i),p,'k-',SAHR_3(:,i),p,'r-',SALR_brut(:,i),p,'b-',SALR_3(:,i),p,'g-')
        xlabel('SA (g.kg^{-1})')
        if sum(~isnan(SAHR_3(:,i)))
            xlim([min(SAHR_3(:,i))-0.1 max(SAHR_3(:,i))+0.1])
        else
            xlim([33 35])
        end
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,3)
        plot(DensHR_brut(:,i),p,'k-',DensHR_3(:,i),p,'r-',DensLR_brut(:,i),p,'b-',DensLR_3(:,i),p,'g-')
        xlabel('\rho_{0} (kg.m^{-3})')
        if sum(~isnan(DensHR_3(:,i)))
            xlim([min(DensHR_3(:,i))-0.1 max(DensHR_3(:,i))+0.1])
        else
        xlim([26 27.5])
        end
        legend('HR RAW','HR GF','LR RAW','LR GF','Location','southwest')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,4)
        plot(CTHR_brut(:,i)-CTHR_3(:,i),p,'b-',CTLR_brut(:,i)-CTLR_3(:,i),p,'r-')
        T=find(~isnan(CTHR_brut(:,i)) & ~isnan(CTHR_3(:,i)));
        ecart_T(i,1)=rms(CTHR_brut(T,i)-CTHR_3(T,i));
        T=find(~isnan(CTLR_brut(:,i)) & ~isnan(CTLR_3(:,i)));
        ecart_L(i,1)=rms(CTLR_brut(T,i)-CTLR_3(T,i));
        anno=['rms HR: ',num2str(ecart_T(i,1),'%.4f') ' rms LR: ' num2str(ecart_L(i,1),'%.4f') ];
        title(anno);
        xlabel('CT RAW - CT GF')
        ylabel('Depth (m)')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,5)
        plot(SAHR_brut(:,i)-SAHR_3(:,i),p,'b-',SALR_brut(:,i)-SALR_3(:,i),p,'r-')
        T=find(~isnan(SAHR_brut(:,i)) & ~isnan(SAHR_3(:,i)));
        ecart_T(i,1)=rms(SAHR_brut(T,i)-SAHR_3(T,i));
        T=find(~isnan(SALR_brut(:,i)) & ~isnan(SALR_3(:,i)));
        ecart_L(i,1)=rms(SALR_brut(T,i)-SALR_3(T,i));
        anno=['rms HR: ',num2str(ecart_T(i,1),'%.4f') ' rms LR: ' num2str(ecart_L(i,1),'%.4f') ];
        title(anno);
        xlabel('SA RAW - SA GF')
        set(gca, 'YDir', 'reverse')
        grid on
        
        subplot(2,3,6)
        plot(DensHR_brut(:,i)-DensHR_3(:,i),p,'b-',DensLR_brut(:,i)-DensLR_3(:,i),p,'r-')
        T=find(~isnan(DensHR_brut(:,i)) & ~isnan(DensHR_3(:,i)));
        ecart_T(i,1)=rms(DensHR_brut(T,i)-DensHR_3(T,i));
        T=find(~isnan(DensLR_brut(:,i)) & ~isnan(DensLR_3(:,i)));
        ecart_L(i,1)=rms(DensLR_brut(T,i)-DensLR_3(T,i));
        anno=['rms HR: ',num2str(ecart_T(i,1),'%.4f') ' rms LR: ' num2str(ecart_L(i,1),'%.4f') ];
        title(anno);
        xlabel('Dens RAW - Dens GF')
        set(gca, 'YDir', 'reverse')
        grid on
        legend('HR RAW - HR GF','LR RAW - LR GF','Location','southeast')
        
        suptitle(['HR & LR RAW VS STEP 3 (GF) - ' nc_att.smru_platform_code ' - Profile # '  num2str(i) ' - # points LR: ' num2str(nbpt_brut(i))]);
        set(gcf,'PaperPositionMode','auto')
        print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/HR_LR/step3/' nc_att.smru_platform_code '_profil_' num2str(i)]);
    end
end
%% comparaison des profils basse resolution en fonction du nombre de points par profils
% CTLR_opt = csvread('CT_opt_interp.csv');
% DensLR_opt = csvread('DensLR_raw_33_interp.txt');
mkdir([dirplot nc_att.smru_platform_code '/comparison/mean_std/']);
clf
p=Pres(:,1);
% on classe les profils en fonction du nombre de points
% seuil à 10 points
I=find(nbpt_brut<10);
J=find(nbpt_brut>=10);
diff_temp=CTHR_brut-CTLR_brut;
%diff_temp_opt = CTHR_brut-CTLR_opt;
diff_sal=SAHR_brut-SALR_brut;
diff_dens=DensHR_brut-DensLR_brut;
%diff_dens_opt = DensHR_brut-DensLR_opt;
% 
% subplot(1,2,1)
% plot(nanmean(diff_dens'),p,nanmean(diff_dens(:,I)'),p,nanmean(diff_dens(:,J)'),p)
% hold on 
% plot(nanmean(diff_dens_opt'),p,nanmean(diff_dens_opt(:,I)'),p,nanmean(diff_dens_opt(:,J)'),p)
% hold off
% grid on
% xlabel('Mean diff \rho_{0}')
% ylabel('Depth (m)')
% xlim([-0.01 0.01])
% set(gca, 'YDir', 'reverse')
% legend('linear LR all','linear LR < 10pts','linear LR >= 10pts','opt LR all','opt LR < 10pts','opt LR >= 10pts','location','southwest')
% 
% subplot(1,2,2)
% plot(nanstd(diff_dens'),p,nanstd(diff_dens(:,I)'),p,nanstd(diff_dens(:,J)'),p)
% hold on 
% plot(nanstd(diff_dens_opt'),p,nanstd(diff_dens_opt(:,I)'),p,nanstd(diff_dens_opt(:,J)'),p)
% hold off
% grid on
% xlabel('STD diff \rho_{0}')
% ylabel('Depth (m)')
% xlim([0 0.04])
% set(gca, 'YDir', 'reverse')

subplot(2,3,1)
plot(nanmean(diff_temp(10:end,:)'),p(10:end),'-g',nanmean(diff_temp(10:end,I)'),p(10:end),'-k')
xct = xlim;

subplot(2,3,2)
plot(nanmean(diff_sal(10:end,:)'),p(10:end),'-g',nanmean(diff_sal(10:end,I)'),p(10:end),'-k')
xsa = xlim;

subplot(2,3,3)
plot(nanmean(diff_dens(10:end,:)'),p(10:end),'-g',nanmean(diff_dens(10:end,I)'),p(10:end),'-k')
xdens = xlim;

subplot(2,3,4)
plot(nanstd(diff_temp(10:end,:)'),p(10:end),'-g',nanstd(diff_temp(10:end,I)'),p(10:end),'-k')
xstdct=xlim;

subplot(2,3,5)
plot(nanstd(diff_sal(10:end,:)'),p(10:end),'-g',nanstd(diff_sal(10:end,I)'),p(10:end),'-k')
xstdsa= xlim;

subplot(2,3,6)
plot(nanstd(diff_dens(10:end,:)'),p(10:end),'-g',nanstd(diff_dens(10:end,I)'),p(10:end),'-k')
xstddens = xlim;


% set axis limits 
xct = [-max(abs(xct)) max(abs(xct))];
xsa = [-max(abs(xsa)) max(abs(xsa))];
xdens = [-max(abs(xdens)) max(abs(xdens))];
clf;
subplot(2,3,1)
plot(nanmean(diff_temp'),p,'-k',nanmean(diff_temp(:,I)'),p,'-b',nanmean(diff_temp(:,J)'),p,'-r')
grid on
xlabel('Mean diff CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
legend({'LR all','LR < 10pts','LR >= 10pts'},'location','southwest','fontsize',5)
xlim(xct)

subplot(2,3,2)
plot(nanmean(diff_sal'),p,'-k',nanmean(diff_sal(:,I)'),p,'-b',nanmean(diff_sal(:,J)'),p,'-r')
xlabel('Mean diff SA (g.kg^{-1})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xsa)

subplot(2,3,3)
plot(nanmean(diff_dens'),p,'-k',nanmean(diff_dens(:,I)'),p,'-b',nanmean(diff_dens(:,J)'),p,'-r')
xlabel('Mean diff \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xdens)

subplot(2,3,4)
plot(nanstd(diff_temp'),p,'-k',nanstd(diff_temp(:,I)'),p,'-b',nanstd(diff_temp(:,J)'),p,'-r')
xlabel('STD diff CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
grid on
xlim(xstdct)

subplot(2,3,5)
plot(nanstd(diff_sal'),p,'-k',nanstd(diff_sal(:,I)'),p,'-b',nanstd(diff_sal(:,J)'),p,'-r')
xlabel('STD diff SA (g.kg^{-1})') 
set(gca, 'YDir', 'reverse')
xlim(xstdsa)
grid on

subplot(2,3,6)
plot(nanstd(diff_dens'),p,'-k',nanstd(diff_dens(:,I)'),p,'-b',nanstd(diff_dens(:,J)'),p,'-r')
xlabel('STD diff \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xstddens)

dim = [.27 .628 .2 .3];
str = ['LR tot: ' num2str(length(d_date)) ' profiles, LR < 10pts = ' num2str(length(I)/length(d_date)*100,'%.0f') '%, LR >= 10pts = ' num2str(length(J)/length(d_date)*100,'%.0f') '%'];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
suptitle(['individu : ' nc_att.smru_platform_code ' - mean and std HR_{raw} - LR_{raw} with LR threshold = 10 pts']);
% set(gcf,'PaperPositionMode','auto')
print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/mean_std' '/mean_std_nbr_pt_raw_' nc_att.smru_platform_code]);

%% step 1
clf
p=Pres(:,1);
% on classe les profils en fonction du nombre de points
% seuil à 10 points
I=find(nbpt_brut<10);
J=find(nbpt_brut>=10);
diff_temp=CTHR_cor-CTLR_cor;
diff_sal=SAHR_cor-SALR_cor;
diff_dens=DensHR_cor-DensLR_cor;

subplot(2,3,1)
plot(nanmean(diff_temp'),p,'-k',nanmean(diff_temp(:,I)'),p,'-b',nanmean(diff_temp(:,J)'),p,'-r')
grid on
xlabel('Mean diff CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
xlim(xct)
legend({'LR all','LR < 10pts','LR >= 10pts'},'location','southwest','fontsize',5)

subplot(2,3,2)
plot(nanmean(diff_sal'),p,'-k',nanmean(diff_sal(:,I)'),p,'-b',nanmean(diff_sal(:,J)'),p,'-r')
xlabel('Mean diff SA (g.kg^{-1})')
xlim([-0.02 0.02])
set(gca, 'YDir', 'reverse')
grid on
xlim(xsa)

subplot(2,3,3)
plot(nanmean(diff_dens'),p,'-k',nanmean(diff_dens(:,I)'),p,'-b',nanmean(diff_dens(:,J)'),p,'-r')
xlabel('Mean diff \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xdens)

subplot(2,3,4)
plot(nanstd(diff_temp'),p,'-k',nanstd(diff_temp(:,I)'),p,'-b',nanstd(diff_temp(:,J)'),p,'-r')
xlabel('STD diff CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
grid on
xlim(xstdct)

subplot(2,3,5)
plot(nanstd(diff_sal'),p,'-k',nanstd(diff_sal(:,I)'),p,'-b',nanstd(diff_sal(:,J)'),p,'-r')
xlabel('STD diff SA (g.kg^{-1})') 
set(gca, 'YDir', 'reverse')
grid on
xlim(xstdsa)

subplot(2,3,6)
plot(nanstd(diff_dens'),p,'-k',nanstd(diff_dens(:,I)'),p,'-b',nanstd(diff_dens(:,J)'),p,'-r')
xlabel('STD diff \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xstddens)

dim = [.27 .628 .2 .3];
str = ['LR tot: ' num2str(length(d_date)) ' profiles, LR < 10pts = ' num2str(length(I)/length(d_date)*100,'%.0f') '%, LR >= 10pts = ' num2str(length(J)/length(d_date)*100,'%.0f') '%'];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
suptitle(['individu : ' nc_att.smru_platform_code ' - mean and std HR_{TME} - LR_{TME} with LR threshold = 10 pts']);
% set(gcf,'PaperPositionMode','auto')
print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/comparison//mean_std' '/mean_std_nbr_pt_tme_' nc_att.smru_platform_code]);
%% step 2
clf
p=Pres(:,1);
% on classe les profils en fonction du nombre de points
% seuil à 10 points
I=find(nbpt_brut<10);
J=find(nbpt_brut>=10);
diff_temp=CTHR_2-CTLR_2;
diff_sal=SAHR_2-SALR_2;
diff_dens=DensHR_2-DensLR_2;

subplot(2,3,1)
plot(nanmean(diff_temp'),p,'-k',nanmean(diff_temp(:,I)'),p,'-b',nanmean(diff_temp(:,J)'),p,'-r')
grid on
xlabel('Mean diff CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
xlim(xct)
legend({'LR all','LR < 10pts','LR >= 10pts'},'location','southwest','fontsize',5)

subplot(2,3,2)
plot(nanmean(diff_sal'),p,'-k',nanmean(diff_sal(:,I)'),p,'-b',nanmean(diff_sal(:,J)'),p,'-r')
xlabel('Mean diff SA (g.kg^{-1})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xsa)

subplot(2,3,3)
plot(nanmean(diff_dens'),p,'-k',nanmean(diff_dens(:,I)'),p,'-b',nanmean(diff_dens(:,J)'),p,'-r')
xlabel('Mean diff \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xdens)

subplot(2,3,4)
plot(nanstd(diff_temp'),p,'-k',nanstd(diff_temp(:,I)'),p,'-b',nanstd(diff_temp(:,J)'),p,'-r')
xlabel('STD diff CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
grid on
xlim(xstdct)

subplot(2,3,5)
plot(nanstd(diff_sal'),p,'-k',nanstd(diff_sal(:,I)'),p,'-b',nanstd(diff_sal(:,J)'),p,'-r')
xlabel('STD diff SA (g.kg^{-1})') 
set(gca, 'YDir', 'reverse')
grid on
xlim(xstdsa)

subplot(2,3,6)
plot(nanstd(diff_dens'),p,'-k',nanstd(diff_dens(:,I)'),p,'-b',nanstd(diff_dens(:,J)'),p,'-r')
xlabel('STD diff \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xstddens)


dim = [.27 .628 .2 .3];
str = ['LR tot: ' num2str(length(d_date)) ' profiles, LR < 10pts = ' num2str(length(I)/length(d_date)*100,'%.0f') '%, LR >= 10pts = ' num2str(length(J)/length(d_date)*100,'%.0f') '%'];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
suptitle(['individu : ' nc_att.smru_platform_code ' - mean and std HR_{DIR} - LR_{DIR} with LR threshold = 10 pts']);
% set(gcf,'PaperPositionMode','auto')
print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/comparison//mean_std' '/mean_std_nbr_pt_dir_' nc_att.smru_platform_code] );
%% step 3
clf
p=Pres(:,1);
% on classe les profils en fonction du nombre de points
% seuil à 10 points
I=find(nbpt_brut<10);
J=find(nbpt_brut>=10);
diff_temp=CTHR_3-CTLR_3;
diff_sal=SAHR_3-SALR_3;
diff_dens=DensHR_3-DensLR_3;

subplot(2,3,1)
plot(nanmean(diff_temp'),p,'-k',nanmean(diff_temp(:,I)'),p,'-b',nanmean(diff_temp(:,J)'),p,'-r')
grid on
xlabel('Mean diff CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
xlim(xct)
legend({'LR all','LR < 10pts','LR >= 10pts'},'location','southwest','fontsize',5)

subplot(2,3,2)
plot(nanmean(diff_sal'),p,'-k',nanmean(diff_sal(:,I)'),p,'-b',nanmean(diff_sal(:,J)'),p,'-r')
xlabel('Mean diff SA (g.kg^{-1})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xsa)

subplot(2,3,3)
plot(nanmean(diff_dens'),p,'-k',nanmean(diff_dens(:,I)'),p,'-b',nanmean(diff_dens(:,J)'),p,'-r')
xlabel('Mean diff \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xdens)

subplot(2,3,4)
plot(nanstd(diff_temp'),p,'-k',nanstd(diff_temp(:,I)'),p,'-b',nanstd(diff_temp(:,J)'),p,'-r')
xlabel('STD diff CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
grid on
xlim(xstdct)

subplot(2,3,5)
plot(nanstd(diff_sal'),p,'-k',nanstd(diff_sal(:,I)'),p,'-b',nanstd(diff_sal(:,J)'),p,'-r')
xlabel('STD diff SA (g.kg^{-1})') 
set(gca, 'YDir', 'reverse')
grid on
xlim(xstdsa)

subplot(2,3,6)
plot(nanstd(diff_dens'),p,'-k',nanstd(diff_dens(:,I)'),p,'-b',nanstd(diff_dens(:,J)'),p,'-r')
xlabel('STD diff \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xstddens)


dim = [.27 .628 .2 .3];
str = ['LR tot: ' num2str(length(d_date)) ' profiles, LR < 10pts = ' num2str(length(I)/length(d_date)*100,'%.0f') '%, LR >= 10pts = ' num2str(length(J)/length(d_date)*100,'%.0f') '%'];
annotation('textbox',dim,'String',str,'FitBoxToText','on');
suptitle(['individu : ' nc_att.smru_platform_code ' - mean and std HR_{GF} - LR_{GF} with LR threshold = 10 pts']);
% set(gcf,'PaperPositionMode','auto')
% print('-dpng',[dirplot nc_att.smru_platform_code '/comparison//mean_std' '/mean_std_nbr_pt_gf_' nc_att.smru_platform_code],'-r400');
print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/comparison//mean_std' '/mean_std_nbr_pt_gf_' nc_att.smru_platform_code]);
%% comparaison des profils basse resolution en fonction du nombre de points par profils RAW vs STEP3

clf
% on classe les profils en fonction du nombre de points
% seuil à 10 points
I=find(nbpt_brut<10);
J=find(nbpt_brut>=10);
p=Pres(:,1);

diff_temp=CTLR_brut-CTLR_3;
diff_sal=SALR_brut-SALR_3;
diff_dens=DensLR_brut-DensLR_3;

diff_tempHR=CTHR_brut-CTHR_3;
diff_salHR=SAHR_brut-SAHR_3;
diff_densHR=DensHR_brut-DensHR_3;


subplot(2,3,1)
plot(nanmean(diff_tempHR(10:end,:)'),p(10:end),'-g')
xct = xlim;

subplot(2,3,2)
plot(nanmean(diff_salHR(10:end,:)'),p(10:end),'-g')
xsa = xlim;

subplot(2,3,3)
plot(nanmean(diff_densHR(10:end,:)'),p(10:end),'-g')
xdens = xlim;

subplot(2,3,4)
plot(nanstd(diff_tempHR(10:end,:)'),p(10:end),'-g')
xstdct=xlim;

subplot(2,3,5)
plot(nanstd(diff_salHR(10:end,:)'),p(10:end),'-g')
xstdsa= xlim;

subplot(2,3,6)
plot(nanstd(diff_densHR(10:end,:)'),p(10:end),'-g')
xstddens = xlim;


% set axis limits 
xct = [-max(abs(xct)) max(abs(xct))];
xsa = [-max(abs(xsa)) max(abs(xsa))];
xdens = [-max(abs(xdens)) max(abs(xdens))];

% plot with set axis lim
subplot(2,3,1)
plot(nanmean(diff_tempHR'),p,'-g',nanmean(diff_temp'),p,'-k',nanmean(diff_temp(:,I)'),p,'-b',nanmean(diff_temp(:,J)'),p,'-r')
grid on
xlabel('Mean diff CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
legend({'HR','LR all','LR < 10pts','LR >= 10pts'},'location','southwest','fontsize',5)
xlim(xct)

subplot(2,3,2)
plot(nanmean(diff_salHR'),p,'-g',nanmean(diff_sal'),p,'-k',nanmean(diff_sal(:,I)'),p,'-b',nanmean(diff_sal(:,J)'),p,'-r')
xlabel('Mean diff SA (g.kg^{-1})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xsa)

subplot(2,3,3)
plot(nanmean(diff_densHR'),p,'-g',nanmean(diff_dens'),p,'-k',nanmean(diff_dens(:,I)'),p,'-b',nanmean(diff_dens(:,J)'),p,'-r')
xlabel('Mean diff \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xdens)

subplot(2,3,4)
plot(nanstd(diff_tempHR'),p,'-g',nanstd(diff_temp'),p,'-k',nanstd(diff_temp(:,I)'),p,'-b',nanstd(diff_temp(:,J)'),p,'-r')
xlabel('STD diff CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
grid on
xlim(xstdct)

subplot(2,3,5)
plot(nanstd(diff_salHR'),p,'-g',nanstd(diff_sal'),p,'-k',nanstd(diff_sal(:,I)'),p,'-b',nanstd(diff_sal(:,J)'),p,'-r')
xlabel('STD diff SA (g.kg^{-1})') 
set(gca, 'YDir', 'reverse')
grid on
xlim(xstdsa)

subplot(2,3,6)
plot(nanstd(diff_densHR'),p,'-g',nanstd(diff_dens'),p,'-k',nanstd(diff_dens(:,I)'),p,'-b',nanstd(diff_dens(:,J)'),p,'-r')
xlabel('STD diff \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xstddens)

dim = [.27 .628 .2 .3];
str = ['LR tot: ' num2str(length(d_date)) ' profiles, LR < 10pts = ' num2str(length(I)/length(d_date)*100,'%.0f') '%, LR >= 10pts = ' num2str(length(J)/length(d_date)*100,'%.0f') '%'];
annotation('textbox',dim,'String',str,'FitBoxToText','on');

suptitle([ nc_att.smru_platform_code ' - mean and std difference RAW - STEP 3 (GF)']);

%set(gcf,'PaperPositionMode','auto')
print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/mean_std/' 'mean_std_diff_raw_GF_' nc_att.smru_platform_code] );

%% comparaison des profils basse resolution en fonction du nombre de points par profils RAW vs STEP2

clf
% on classe les profils en fonction du nombre de points
% seuil à 10 points
I=find(nbpt_brut<10);
J=find(nbpt_brut>=10); 
p=Pres(:,1);

diff_temp=CTLR_brut-CTLR_2;
diff_sal=SALR_brut-SALR_2;
diff_dens=DensLR_brut-DensLR_2;

diff_tempHR=CTHR_brut-CTHR_2;
diff_salHR=SAHR_brut-SAHR_2;
diff_densHR=DensHR_brut-DensHR_2;


subplot(2,3,1)
plot(nanmean(diff_tempHR'),p,'-g',nanmean(diff_temp'),p,'-k',nanmean(diff_temp(:,I)'),p,'-b',nanmean(diff_temp(:,J)'),p,'-r')
grid on
xlabel('Mean diff CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
legend({'HR','LR all','LR < 10pts','LR >= 10pts'},'location','southwest','fontsize',5)
xlim(xct)

subplot(2,3,2)
plot(nanmean(diff_salHR'),p,'-g',nanmean(diff_sal'),p,'-k',nanmean(diff_sal(:,I)'),p,'-b',nanmean(diff_sal(:,J)'),p,'-r')
xlabel('Mean diff SA (g.kg^{-1})')
xlim(xsa)
set(gca, 'YDir', 'reverse')
grid on

subplot(2,3,3)
plot(nanmean(diff_densHR'),p,'-g',nanmean(diff_dens'),p,'-k',nanmean(diff_dens(:,I)'),p,'-b',nanmean(diff_dens(:,J)'),p,'-r')
xlabel('Mean diff \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xdens)

subplot(2,3,4)
plot(nanstd(diff_tempHR'),p,'-g',nanstd(diff_temp'),p,'-k',nanstd(diff_temp(:,I)'),p,'-b',nanstd(diff_temp(:,J)'),p,'-r')
xlabel('STD diff CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
grid on
xlim(xstdct)

subplot(2,3,5)
plot(nanstd(diff_salHR'),p,'-g',nanstd(diff_sal'),p,'-k',nanstd(diff_sal(:,I)'),p,'-b',nanstd(diff_sal(:,J)'),p,'-r')
xlabel('STD diff SA (g.kg^{-1})') 
set(gca, 'YDir', 'reverse')
xlim(xstdsa)
grid on

subplot(2,3,6)
plot(nanstd(diff_densHR'),p,'-g',nanstd(diff_dens'),p,'-k',nanstd(diff_dens(:,I)'),p,'-b',nanstd(diff_dens(:,J)'),p,'-r')
xlabel('STD diff \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xstddens)
 
dim = [.27 .628 .2 .3];
str = ['LR tot: ' num2str(length(d_date)) ' profiles, LR < 10pts = ' num2str(length(I)/length(d_date)*100,'%.0f') '%, LR >= 10pts = ' num2str(length(J)/length(d_date)*100,'%.0f') '%'];
annotation('textbox',dim,'String',str,'FitBoxToText','on');

suptitle([ nc_att.smru_platform_code ' - mean and std difference RAW - STEP 2 (DIR)']);

%set(gcf,'PaperPositionMode','auto')
print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/mean_std/' 'mean_std_diff_raw_DIR_' nc_att.smru_platform_code] );

%% comparaison des profils basse resolution en fonction du nombre de points par profils RAW vs STEP1

clf
% on classe les profils en fonction du nombre de points
% seuil à 10 points
I=find(nbpt_brut<10);
J=find(nbpt_brut>=10);
p=Pres(:,1);

diff_temp=CTLR_brut-CTLR_cor;
diff_sal=SALR_brut-SALR_cor;
diff_dens=DensLR_brut-DensLR_cor;

diff_tempHR=CTHR_brut-CTHR_cor;
diff_salHR=SAHR_brut-SAHR_cor;
diff_densHR=DensHR_brut-DensHR_cor;

subplot(2,3,1)
plot(nanmean(diff_tempHR'),p,'-g',nanmean(diff_temp'),p,'-k',nanmean(diff_temp(:,I)'),p,'-b',nanmean(diff_temp(:,J)'),p,'-r')
grid on
xlabel('Mean diff CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
legend({'HR','LR all','LR < 10pts','LR >= 10pts'},'location','southwest','fontsize',5)
xlim(xct)

subplot(2,3,2)
plot(nanmean(diff_salHR'),p,'-g',nanmean(diff_sal'),p,'-k',nanmean(diff_sal(:,I)'),p,'-b',nanmean(diff_sal(:,J)'),p,'-r')
xlabel('Mean diff SA (g.kg^{-1})')
xlim(xsa)
set(gca, 'YDir', 'reverse')
grid on

subplot(2,3,3)
plot(nanmean(diff_densHR'),p,'-g',nanmean(diff_dens'),p,'-k',nanmean(diff_dens(:,I)'),p,'-b',nanmean(diff_dens(:,J)'),p,'-r')
xlabel('Mean diff \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xdens)

subplot(2,3,4)
plot(nanstd(diff_tempHR'),p,'-g',nanstd(diff_temp'),p,'-k',nanstd(diff_temp(:,I)'),p,'-b',nanstd(diff_temp(:,J)'),p,'-r')
xlabel('STD diff CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
grid on
xlim(xstdct)

subplot(2,3,5)
plot(nanstd(diff_salHR'),p,'-g',nanstd(diff_sal'),p,'-k',nanstd(diff_sal(:,I)'),p,'-b',nanstd(diff_sal(:,J)'),p,'-r')
xlabel('STD diff SA (g.kg^{-1})') 
set(gca, 'YDir', 'reverse')
xlim(xstdsa)
grid on

subplot(2,3,6)
plot(nanstd(diff_densHR'),p,'-g',nanstd(diff_dens'),p,'-k',nanstd(diff_dens(:,I)'),p,'-b',nanstd(diff_dens(:,J)'),p,'-r')
xlabel('STD diff \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xstddens)

dim = [.27 .628 .2 .3];
str = ['LR tot: ' num2str(length(d_date)) ' profiles, LR < 10pts = ' num2str(length(I)/length(d_date)*100,'%.0f') '%, LR >= 10pts = ' num2str(length(J)/length(d_date)*100,'%.0f') '%'];
annotation('textbox',dim,'String',str,'FitBoxToText','on');

suptitle([ nc_att.smru_platform_code ' - mean and std difference RAW - STEP 1 (TME)']);


print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/mean_std/' 'mean_std_diff_raw_TME_' nc_att.smru_platform_code] );

%% transects

colormap default
    Tcontour=[-1:2:15];
    Scontour=[33.6:0.5:34.6];
    Dcontour=[26.5:0.5:27.5];
    pmax=500;
    [Xi,Yi]=meshgrid(d_date-d_date(1),p(1:pmax));
    
    % figure CT HR raw
    clf
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,CTHR_brut(1:pmax,:))
    shading('interp')
    contour(Xi,Yi,CTHR_brut(1:pmax,:),Tcontour,'linecolor','k')
    xct = caxis;
    c=colorbar('northoutside');
    ylabel(c,'CT (°C)')
    xlabel('days')
    ylabel('Depth (m)')
    title([ nc_att.smru_platform_code ' - CT HR RAW (°C)'], 'fontsize', 18);
    grid on
    set(gca, 'YDir', 'reverse')
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/' nc_att.smru_platform_code '_transec_tempHR_raw']);

    % figure SA HR raw
    clf
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,SAHR_brut(1:pmax,:))
    shading('interp')
    contour(Xi,Yi,SAHR_brut(1:pmax,:),Scontour,'linecolor','k')
    xsa=caxis;
    c=colorbar('northoutside');
    ylabel(c,'SA (g.kg^{-1})')
    xlabel('days')
    ylabel('Depth (m)')
    title([ nc_att.smru_platform_code ' - SA HR RAW (g.kg^{-1})'], 'fontsize', 18);
    grid on
    set(gca, 'YDir', 'reverse')
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/' nc_att.smru_platform_code '_transec_salHR_raw']);
    
    % figure CT HR 3
    clf
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,CTHR_3(1:pmax,:))
    shading('interp')
    contour(Xi,Yi,CTHR_3(1:pmax,:),Tcontour,'linecolor','k')
    caxis(xct)
    c=colorbar('northoutside');
    ylabel(c,'CT (°C)')
    xlabel('days')
    ylabel('Depth (m)')
    title([ nc_att.smru_platform_code ' - CT HR 3 (°C)'], 'fontsize', 18);
    grid on
    set(gca, 'YDir', 'reverse')
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/' nc_att.smru_platform_code '_transec_tempHR_3']);

    % figure SA HR 3
    clf
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,SAHR_3(1:pmax,:))
    shading('interp')
    contour(Xi,Yi,SAHR_3(1:pmax,:),Scontour,'linecolor','k')
    caxis(xsa)
    c=colorbar('northoutside');
    ylabel(c,'SA (g.kg^{-1})')
    xlabel('days')
    ylabel('Depth (m)')
    title([ nc_att.smru_platform_code ' - SA HR 3 (g.kg^{-1})'], 'fontsize', 18);
    grid on
    set(gca, 'YDir', 'reverse')
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/' nc_att.smru_platform_code '_transec_salHR_3']);
    
    % figure CT HR raw
    clf
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,CTLR_brut(1:pmax,:))
    shading('interp')
    contour(Xi,Yi,CTLR_brut(1:pmax,:),Tcontour,'linecolor','k')
    caxis(xct)
    c=colorbar('northoutside');
    ylabel(c,'CT (°C)')
    xlabel('days')
    ylabel('Depth (m)')
    title([ nc_att.smru_platform_code ' - CT LR RAW (°C)'], 'fontsize', 18);
    grid on
    set(gca, 'YDir', 'reverse')
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/' nc_att.smru_platform_code '_transec_tempLR_raw']);

    % figure SA HR raw
    clf
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,SALR_brut(1:pmax,:))
    shading('interp')
    contour(Xi,Yi,SALR_brut(1:pmax,:),Scontour,'linecolor','k')
    caxis(xsa)
    c=colorbar('northoutside');
    ylabel(c,'SA (g.kg^{-1})')
    xlabel('days')
    ylabel('Depth (m)')
    title([ nc_att.smru_platform_code ' - SA LR RAW (g.kg^{-1})'], 'fontsize', 18);
    grid on
    set(gca, 'YDir', 'reverse')
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/' nc_att.smru_platform_code '_transec_salLR_raw']);
    
    % figure CT LR 3
    clf
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,CTLR_3(1:pmax,:))
    shading('interp')
    contour(Xi,Yi,CTLR_3(1:pmax,:),Tcontour,'linecolor','k')
    caxis(xct)
    c=colorbar('northoutside');
    ylabel(c,'CT (°C)')
    xlabel('days')
    ylabel('Depth (m)')
    title([ nc_att.smru_platform_code ' - CT LR 3 (°C)'], 'fontsize', 18);
    grid on
    set(gca, 'YDir', 'reverse')
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/' nc_att.smru_platform_code '_transec_tempLR_3']);

    % figure SA HR 3
    clf
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,SALR_3(1:pmax,:))
    shading('interp')
    contour(Xi,Yi,SALR_3(1:pmax,:),Scontour,'linecolor','k')
    caxis(xsa)
    c=colorbar('northoutside');
    ylabel(c,'SA (g.kg^{-1})')
    xlabel('days')
    ylabel('Depth (m)')
    title([ nc_att.smru_platform_code ' - SA LR 3 (g.kg^{-1})'], 'fontsize', 18);
    grid on
    set(gca, 'YDir', 'reverse')
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/' nc_att.smru_platform_code '_transec_salLR_3']);
    
    % figure Dens HR raw
    clf
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,DensHR_brut(1:pmax,:))
    shading('interp')
    contour(Xi,Yi,DensHR_brut(1:pmax,:),Dcontour,'linecolor','k')
    xdens = caxis;
    c=colorbar('northoutside');
    ylabel(c,' \rho_{0}  (kg.m^{-3})')
    xlabel('days')
    ylabel('Depth (m)')
    title([ nc_att.smru_platform_code ' - Dens HR raw (kg.m^{-3})'], 'fontsize', 18);
    grid on
    set(gca, 'YDir', 'reverse')
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/' nc_att.smru_platform_code '_transec_densHR_raw']);
    
    % figure Dens LR raw
    clf
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,DensLR_brut(1:pmax,:))
    shading('interp')
    contour(Xi,Yi,DensLR_brut(1:pmax,:),Dcontour,'linecolor','k')
    caxis(xdens)
    c=colorbar('northoutside');
    ylabel(c,' \rho_{0}  (kg.m^{-3})')
    xlabel('days')
    ylabel('Depth (m)')
    title([ nc_att.smru_platform_code ' - Dens LR raw (kg.m^{-3})'], 'fontsize', 18);
    grid on
    set(gca, 'YDir', 'reverse')
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/' nc_att.smru_platform_code '_transec_densLR_raw']);
    
    % figure Dens LR 3
    clf
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,DensLR_3(1:pmax,:))
    shading('interp')
    contour(Xi,Yi,DensLR_3(1:pmax,:),Dcontour,'linecolor','k')
    caxis(xdens)
    c=colorbar('northoutside');
    ylabel(c,' \rho_{0}  (kg.m^{-3})')
    xlabel('days')
    ylabel('Depth (m)')
    title([ nc_att.smru_platform_code ' - Dens LR 3 (kg.m^{-3})'], 'fontsize', 18);
    grid on
    set(gca, 'YDir', 'reverse')
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/' nc_att.smru_platform_code '_transec_densLR_3']);
    
   
    % figure Dens HR 3
    clf
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,DensHR_3(1:pmax,:))
    shading('interp')
    contour(Xi,Yi,DensHR_3(1:pmax,:),Dcontour,'linecolor','k')
    caxis(xdens)
    c=colorbar('northoutside');
    ylabel(c,' \rho_{0}  (kg.m^{-3})')
    xlabel('days')
    ylabel('Depth (m)')
    title([ nc_att.smru_platform_code ' - Dens HR 3 (kg.m^{-3})'], 'fontsize', 18);
    grid on
    set(gca, 'YDir', 'reverse')
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/' nc_att.smru_platform_code '_transec_densHR_3']);
    

%% comparaison HR cor/ LR cor

    clf
    [Xi,Yi]=meshgrid(d_date-d_date(1),p);

    % compare TEMP LR RAW/ TEMP HR RAW
    diff_temp=CTHR_brut-CTLR_brut;
    pal=parula;
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,diff_temp)
    shading('interp')
    colormap(pal)
    caxis([-0.8 0.8])
    c=colorbar;
    ylabel(c,'Diff CT')
    axis([min(Xi(:)) max(Xi(:)) 0 800]);
    xlabel('Days')
    ylabel('Depth (m)')
    set(gca, 'YDir', 'reverse')
    title(['individu : ' nc_att.smru_platform_code ' Difference CT HR raw - CT LR raw']);
    grid on
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf','-r300',[dirplot nc_att.smru_platform_code '/comparison/' '/difference_temperature_hrraw-lrraw_' nc_att.smru_platform_code]);
    
    % compare TEMP LR 3/ TEMP HR 3
    diff_temp=CTHR_3-CTLR_3;
    pal=parula;
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',... 
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,diff_temp)
    shading('interp')
    colormap(pal)
    caxis([-0.8 0.8])
    c=colorbar;
    ylabel(c,'Diff CT')
    axis([min(Xi(:)) max(Xi(:)) 0 800]);
    xlabel('Days')
    ylabel('Depth (m)')
    set(gca, 'YDir', 'reverse')
    title(['individu : ' nc_att.smru_platform_code ' Difference CT HR 3 - CT LR 3']);
    grid on
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf','-r300',[dirplot nc_att.smru_platform_code '/comparison/' '/difference_temperature_hrr3-lrr3_' nc_att.smru_platform_code]);
    
     % compare S LR RAW/ S HR RAW
    diff_sal=SAHR_brut-SALR_brut;
    pal=parula;
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,diff_sal)
    shading('interp')
    colormap(pal)
    caxis([-0.2 0.2])
    c=colorbar;
    ylabel(c,'Diff SA')
    axis([min(Xi(:)) max(Xi(:)) 0 800]);
    xlabel('Days')
    ylabel('Depth (m)')
    set(gca, 'YDir', 'reverse')
    title(['individu : ' nc_att.smru_platform_code ' Difference SA HR raw - SA LR raw']);
    grid on
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf','-r300',[dirplot nc_att.smru_platform_code '/comparison/' '/difference_sal_hrraw-lrraw_' nc_att.smru_platform_code]);
    
    % compare S LR 3/ S HR 3
    diff_sal=SAHR_3-SALR_3;
    pal=parula;
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',... 
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,diff_sal)
    shading('interp')
    colormap(pal)
    caxis([-0.2 0.2])
    c=colorbar;
    ylabel(c,'Diff SA')
    axis([min(Xi(:)) max(Xi(:)) 0 800]);
    xlabel('Days')
    ylabel('Depth (m)')
    set(gca, 'YDir', 'reverse')
    title(['individu : ' nc_att.smru_platform_code ' Difference SA HR 3 - SA LR 3']);
    grid on
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf','-r300',[dirplot nc_att.smru_platform_code '/comparison/' '/difference_sal_hrr3-lrr3_' nc_att.smru_platform_code]);
    
     % compare D LR RAW/ D HR RAW
    diff_dens=DensHR_brut-DensLR_brut;
    pal=parula;
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,diff_dens)
    shading('interp')
    colormap(pal)
    caxis([-0.1 0.1])
    c=colorbar;
    ylabel(c,'Diff dens')
    axis([min(Xi(:)) max(Xi(:)) 0 800]);
    xlabel('Days')
    ylabel('Depth (m)')
    set(gca, 'YDir', 'reverse')
    title(['individu : ' nc_att.smru_platform_code ' Difference Dens HR raw - Dens LR raw']);
    grid on
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf','-r300',[dirplot nc_att.smru_platform_code '/comparison/' '/difference_dens_hrraw-lrraw_' nc_att.smru_platform_code]);
    
    % compare D LR 3/ D HR 3
    diff_dens=DensHR_3-DensLR_3;
    pal=parula;
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',... 
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,Yi,diff_dens)
    shading('interp')
    colormap(pal)
    caxis([-0.1 0.1])
    c=colorbar;
    ylabel(c,'Diff dens')
    axis([min(Xi(:)) max(Xi(:)) 0 800]);
    xlabel('Days')
    ylabel('Depth (m)')
    set(gca, 'YDir', 'reverse')
    title(['individu : ' nc_att.smru_platform_code ' Difference Dens HR 3 - Dens LR 3']);
    grid on
    set(gcf,'PaperPositionMode','auto')
    print('-bestfit','-dpdf','-r300',[dirplot nc_att.smru_platform_code '/comparison/' '/difference_dens_hrr3-lrr3_' nc_att.smru_platform_code]);

%% RMS calc

% Root mean square calculation for the 3 steps 

%step 1 TME LR
T=[];S=[];R=[];
rms_t=[];
rms_s=[];
rms_r=[];
for i=1:size(CTLR_cor,1)
T=find(~isnan(CTLR_cor(i,:)) & ~isnan(CTLR_brut(i,:)));
rms_t(i)=rms(CTLR_cor(i,T)-CTLR_brut(i,T));

S=find(~isnan(SALR_cor(i,:)) & ~isnan(SALR_brut(i,:)));
rms_s(i)=rms(SALR_cor(i,S)-SALR_brut(i,S));

R=find(~isnan(DensLR_cor(i,:)) & ~isnan(DensLR_brut(i,:)));
rms_r(i)=rms(DensLR_cor(i,R)-DensLR_brut(i,R));
end
RmsL_t1=rms_t';
RmsL_s1=rms_s';
RmsL_r1=rms_r';

%step 1 TME HR
T=[];S=[];R=[];
rms_t=[];
rms_s=[];
rms_r=[];
for i=1:size(CTHR_cor,1)
T=find(~isnan(CTHR_cor(i,:)) & ~isnan(CTHR_brut(i,:)));
rms_t(i)=rms(CTHR_cor(i,T)-CTHR_brut(i,T));

S=find(~isnan(SAHR_cor(i,:)) & ~isnan(SAHR_brut(i,:)));
rms_s(i)=rms(SAHR_cor(i,S)-SAHR_brut(i,S));

R=find(~isnan(DensHR_cor(i,:)) & ~isnan(DensHR_brut(i,:)));
rms_r(i)=rms(DensHR_cor(i,R)-DensHR_brut(i,R));
end
Rms_t1=rms_t';
Rms_s1=rms_s';
Rms_r1=rms_r';


%step 2 DIR LR
T=[];S=[];R=[];
rms_t=[];
rms_s=[];
rms_r=[];
for i=1:size(CTLR_2,1)
T=find(~isnan(CTLR_2(i,:)) & ~isnan(CTLR_brut(i,:)));
rms_t(i)=rms(CTLR_2(i,T)-CTLR_brut(i,T));

S=find(~isnan(SALR_2(i,:)) & ~isnan(SALR_brut(i,:)));
rms_s(i)=rms(SALR_2(i,S)-SALR_brut(i,S));

R=find(~isnan(DensLR_2(i,:)) & ~isnan(DensLR_brut(i,:)));
rms_r(i)=rms(DensLR_2(i,R)-DensLR_brut(i,R));
end
RmsL_t2=rms_t';
RmsL_s2=rms_s';
RmsL_r2=rms_r';

%step 2 DIR HR
T=[];S=[];R=[];
rms_t=[];
rms_s=[];
rms_r=[];
for i=1:size(CTHR_2,1)
T=find(~isnan(CTHR_2(i,:)) & ~isnan(CTHR_brut(i,:)));
rms_t(i)=rms(CTHR_2(i,T)-CTHR_brut(i,T));

S=find(~isnan(SAHR_2(i,:)) & ~isnan(SAHR_brut(i,:)));
rms_s(i)=rms(SAHR_2(i,S)-SAHR_brut(i,S));

R=find(~isnan(DensHR_2(i,:)) & ~isnan(DensHR_brut(i,:)));
rms_r(i)=rms(DensHR_2(i,R)-DensHR_brut(i,R));
end
Rms_t2=rms_t';
Rms_s2=rms_s';
Rms_r2=rms_r';

%step 3 gaussian filter LR
T=[];S=[];R=[];
rms_t=[];
rms_s=[];
rms_r=[];

for i=1:size(CTLR_3,1)
T=find(~isnan(CTLR_3(i,:)) & ~isnan(CTLR_brut(i,:)));
rms_t(i)=rms(CTLR_3(i,T)-CTLR_brut(i,T));

S=find(~isnan(SALR_3(i,:)) & ~isnan(SALR_brut(i,:)));
rms_s(i)=rms(SALR_3(i,S)-SALR_brut(i,S));

R=find(~isnan(DensLR_3(i,:)) & ~isnan(DensLR_brut(i,:)));
rms_r(i)=rms(DensLR_3(i,R)-DensLR_brut(i,R));
end
RmsL_t3=rms_t';
RmsL_s3=rms_s';
RmsL_r3=rms_r';

%step 3 gaussian filter HR
T=[];S=[];R=[];
rms_t=[];
rms_s=[];
rms_r=[];

for i=1:size(CTHR_3,1)
T=find(~isnan(CTHR_3(i,:)) & ~isnan(CTHR_brut(i,:)));
rms_t(i)=rms(CTHR_3(i,T)-CTHR_brut(i,T));

S=find(~isnan(SAHR_3(i,:)) & ~isnan(SAHR_brut(i,:)));
rms_s(i)=rms(SAHR_3(i,S)-SAHR_brut(i,S));

R=find(~isnan(DensHR_3(i,:)) & ~isnan(DensHR_brut(i,:)));
rms_r(i)=rms(DensHR_3(i,R)-DensHR_brut(i,R));
end
Rms_t3=rms_t';
Rms_s3=rms_s';
Rms_r3=rms_r';

clear T S R rms_t rms_s rms_r
%% SAVE RMSs
%sigma
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsL_r1.mat'],'RmsL_r1');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsL_r2.mat'],'RmsL_r2');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsL_r3.mat'],'RmsL_r3');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rms_r1.mat'],'Rms_r1');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rms_r2.mat'],'Rms_r2');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rms_r3.mat'],'Rms_r3');
%ct
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsL_t1.mat'],'RmsL_t1');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsL_t2.mat'],'RmsL_t2');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsL_t3.mat'],'RmsL_t3');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rms_t1.mat'],'Rms_t1');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rms_t2.mat'],'Rms_t2');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rms_t3.mat'],'Rms_t3');
%sa
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsL_s1.mat'],'RmsL_s1');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsL_s2.mat'],'RmsL_s2');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsL_s3.mat'],'RmsL_s3');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rms_s1.mat'],'Rms_s1');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rms_s2.mat'],'Rms_s2');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rms_s3.mat'],'Rms_s3');
%% plot RMS

mkdir([dirplot nc_att.smru_platform_code '/comparison/RMS/']);
p=Pres(:,1);
clf
subplot(2,3,1)
plot(Rms_t3(10:end),p(10:end),'-b')
xct=xlim;
yct=ylim;

subplot(2,3,2)
plot(Rms_s3(10:end),p(10:end),'-b')
xsa=xlim;


subplot(2,3,3)
plot(Rms_r3(10:end),p(10:end),'-b')
xdens=xlim;

subplot(2,3,1)
plot(Rms_t1,p,'-g',Rms_t2,p,'-r',Rms_t3,p,'-b')
grid on
xlabel('RMS RAW - STEP X CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
legend('TME','DIR','GF','location','southeast')
xlim(xct)

subplot(2,3,2)
plot(Rms_s1,p,'-g',Rms_s2,p,'-r',Rms_s3,p,'-b')
xlabel('RMS RAW - STEP X SA (g.kg^{-1})')
xlim(xsa)
set(gca, 'YDir', 'reverse')
grid on
title('HR', 'fontsize',15)

subplot(2,3,3)
plot(Rms_r1,p,'-g',Rms_r2,p,'-r',Rms_r3,p,'-b')
xlabel('RMS RAW - STEP X \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xdens)

subplot(2,3,4)
plot(RmsL_t1,p,'-g',RmsL_t2,p,'-r',RmsL_t3,p,'-b')
grid on
xlabel('RMS RAW - STEP X CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
xlim(xct)
grid on
ylim(yct)

subplot(2,3,5)
plot(RmsL_s1,p,'-g',RmsL_s2,p,'-r',RmsL_s3,p,'-b')
grid on
xlabel('RMS RAW - STEP X SA (g.kg^{-1})')
set(gca, 'YDir', 'reverse')
ylim(yct)
xlim(xsa)
title('LR', 'fontsize',15)

subplot(2,3,6)
plot(RmsL_r1,p,'-g',RmsL_r2,p,'-r',RmsL_r3,p,'-b')
xlabel('RMS RAW - STEP X \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xdens)
ylim(yct)


suptitle([ nc_att.smru_platform_code ' - RMS HR AND LR']);

%set(gcf,'PaperPositionMode','auto')
print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/RMS/' '/mean_std_diff_raw-g_' nc_att.smru_platform_code] );

%% RMS LR vs HR at diff step 
% LR/HR RAW
T=[];S=[];R=[];
rms_t=[];
rms_s=[];
rms_r=[];
for i=1:size(CTLR_brut,1)
T=find(~isnan(CTLR_brut(i,:)) & ~isnan(CTHR_brut(i,:)));
rms_t(i)=rms(CTLR_brut(i,T)-CTHR_brut(i,T));

S=find(~isnan(SALR_brut(i,:)) & ~isnan(SAHR_brut(i,:)));
rms_s(i)=rms(SALR_brut(i,S)-SAHR_brut(i,S));

R=find(~isnan(DensLR_brut(i,:)) & ~isnan(DensHR_brut(i,:)));
rms_r(i)=rms(DensLR_brut(i,R)-DensHR_brut(i,R));
end
RmsX_tb=rms_t';
RmsX_sb=rms_s';
RmsX_rb=rms_r';

% LR/HR TME
T=[];S=[];R=[];
rms_t=[];
rms_s=[];
rms_r=[];
for i=1:size(CTLR_cor,1)
T=find(~isnan(CTLR_cor(i,:)) & ~isnan(CTHR_cor(i,:)));
rms_t(i)=rms(CTLR_cor(i,T)-CTHR_cor(i,T));

S=find(~isnan(SALR_cor(i,:)) & ~isnan(SAHR_cor(i,:)));
rms_s(i)=rms(SALR_cor(i,S)-SAHR_cor(i,S));

R=find(~isnan(DensLR_cor(i,:)) & ~isnan(DensHR_cor(i,:)));
rms_r(i)=rms(DensLR_cor(i,R)-DensHR_cor(i,R));
end
RmsX_t1=rms_t';
RmsX_s1=rms_s';
RmsX_r1=rms_r';

% LR/HR DIR
T=[];S=[];R=[];
rms_t=[];
rms_s=[];
rms_r=[];
for i=1:size(CTLR_2,1)
T=find(~isnan(CTLR_2(i,:)) & ~isnan(CTHR_2(i,:)));
rms_t(i)=rms(CTLR_2(i,T)-CTHR_2(i,T));

S=find(~isnan(SALR_2(i,:)) & ~isnan(SAHR_2(i,:)));
rms_s(i)=rms(SALR_2(i,S)-SAHR_2(i,S));

R=find(~isnan(DensLR_2(i,:)) & ~isnan(DensHR_2(i,:)));
rms_r(i)=rms(DensLR_2(i,R)-DensHR_2(i,R));
end
RmsX_t2=rms_t';
RmsX_s2=rms_s';
RmsX_r2=rms_r';

% LR/HR GF
T=[];S=[];R=[];
rms_t=[];
rms_s=[];
rms_r=[];
for i=1:size(CTLR_3,1)
T=find(~isnan(CTLR_3(i,:)) & ~isnan(CTHR_3(i,:)));
rms_t(i)=rms(CTLR_3(i,T)-CTHR_3(i,T));

S=find(~isnan(SALR_3(i,:)) & ~isnan(SAHR_3(i,:)));
rms_s(i)=rms(SALR_3(i,S)-SAHR_3(i,S));

R=find(~isnan(DensLR_3(i,:)) & ~isnan(DensHR_3(i,:)));
rms_r(i)=rms(DensLR_3(i,R)-DensHR_3(i,R));
end
RmsX_t3=rms_t';
RmsX_s3=rms_s';
RmsX_r3=rms_r';
%% SAVE RMSs LR vs HR
%sigma
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsLRHR_rb.mat'],'RmsX_rb');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsLRHR_r1.mat'],'RmsX_r1');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsLRHR_r2.mat'],'RmsX_r2');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsLRHR_r3.mat'],'RmsX_r3');
% ct
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsLRHR_tb.mat'],'RmsX_tb');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsLRHR_t1.mat'],'RmsX_t1');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsLRHR_t2.mat'],'RmsX_t2');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsLRHR_t3.mat'],'RmsX_t3');
% sa
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsLRHR_sb.mat'],'RmsX_sb');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsLRHR_s1.mat'],'RmsX_s1');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsLRHR_s2.mat'],'RmsX_s2');
save(['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/variables/' nc_att.smru_platform_code '/' nc_att.smru_platform_code '_rmsLRHR_s3.mat'],'RmsX_s3');
%% plot RMS LR_HR

clf;

subplot(1,3,1)
plot(RmsX_tb(10:end),p(10:end),'-g')
xct=xlim;

subplot(1,3,2)
plot(RmsX_sb(10:end),p(10:end),'-g')
xsa=xlim;

subplot(1,3,3)
plot(RmsX_rb(10:end),p(10:end),'-g')
xdens=xlim;



subplot(1,3,1)
plot(RmsX_tb,p,'-g',RmsX_t1,p,'-r',RmsX_t2,p,'-b',RmsX_t3,p,'-k')
grid on
xlabel('RMS LR - HR CT (°C)')
ylabel('Depth (m)')
set(gca, 'YDir', 'reverse')
legend('RAW','TME','DIR','GF','location','southeast')
xlim(xct)

subplot(1,3,2)
plot(RmsX_sb,p,'-g',RmsX_s1,p,'-r',RmsX_s2,p,'-b',RmsX_s3,p,'-k')
xlabel('RMS LR - HR SA (g.kg^{-1})')
xlim(xsa)
set(gca, 'YDir', 'reverse')
grid on


subplot(1,3,3)
plot(RmsX_rb,p,'-g',RmsX_r1,p,'-r',RmsX_r2,p,'-b',RmsX_r3,p,'-k')
xlabel('RMS LR - HR \rho_{0} (kg.m^{-3})')
set(gca, 'YDir', 'reverse')
grid on
xlim(xdens)

suptitle([ nc_att.smru_platform_code ' - RMS HR vs LR']);

%set(gcf,'PaperPositionMode','auto')
print('-fillpage','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/RMS/' '/RMS_LR_HR' nc_att.smru_platform_code] );

%% Histogram of maximal dives depth LR +- 10pts and HR

pL_sel = PresL(:,ind_lr);
B = ~isnan(pL_sel);
% indices
ind_pL_max = arrayfun(@(x) find(B(:, x), 1, 'last'), 1:size(pL_sel, 2));
% values
pL_max = arrayfun(@(x,y) pL_sel(x,y), ind_pL_max, 1:size(pL_sel, 2));
 
I=find(nbpt_brut<10);
J=find(nbpt_brut>=10);

clf;
h1=histogram([pL_max(J),pL_max(I)],30,'facealpha',.5);
%pbaspect([1 1 1])
hold on 
h2=histogram(pL_max(J),30,'facealpha',.5);
h3=histogram(pL_max(I),30,'facealpha',.5);
hold off
grid on
hleg=legend([h1,h2,h3],'LR','LR>=10pts','LR<10pts');

title(['BOTTOM PRESSURE (dbar)'])
%set(gcf,'PaperPositionMode','auto')
print('-bestfit','-dpdf',[dirplot nc_att.smru_platform_code '/comparison/' '/depths_' nc_att.smru_platform_code] );

%% Identification of bugged profils : vertical sigma0 rms > threshold between LR and HR

for i=1:size(DensHR_brut,2)
        T=find(~isnan(DensHR_brut(:,i)) & ~isnan(DensLR_brut(:,i)));
        ecart_D(i,1)=rms(DensHR_brut(T,i)-DensLR_brut(T,i));
        ecart_T(i,1)=rms(CTHR_brut(T,i)-CTLR_brut(T,i));
        ecart_S(i,1)=rms(SAHR_brut(T,i)-SALR_brut(T,i));
end
I=find(nbpt_brut<10);
J=find(nbpt_brut>=10);

flag = find(ecart_S>0.06);
x = intersect(flag,I);
y = intersect(flag,J);
p=100.*[0.06/100 length(flag)/length(ecart_S) length(x)/length(I) length(y)/length(J)];

flag = find(ecart_S>0.05);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.05/100 length(flag)/length(ecart_S) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_S>0.04);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.04/100 length(flag)/length(ecart_S) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_S>0.03);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.03/100 length(flag)/length(ecart_S) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_S>0.02);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.02/100 length(flag)/length(ecart_S) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_S>0.01);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.01/100 length(flag)/length(ecart_S) length(x)/length(I) length(y)/length(J)]];


P_sa=array2table(p,'VariableNames',{'threshold_SA','percent_LR_tot','LR_less_10','LR_more_10'})

flag = find(ecart_D>0.06);
x = intersect(flag,I);
y = intersect(flag,J);
p=100.*[0.06/100 length(flag)/length(ecart_D) length(x)/length(I) length(y)/length(J)];

flag = find(ecart_D>0.05);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.05/100 length(flag)/length(ecart_D) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_D>0.04);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.04/100 length(flag)/length(ecart_D) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_D>0.03);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.03/100 length(flag)/length(ecart_D) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_D>0.02);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.02/100 length(flag)/length(ecart_D) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_D>0.01);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.01/100 length(flag)/length(ecart_D) length(x)/length(I) length(y)/length(J)]];

P_dens=array2table(p,'VariableNames',{'threshold_sigma0','percent_LR_tot','LR_less_10','LR_more_10'})

flag = find(ecart_T>0.5);
x = intersect(flag,I);
y = intersect(flag,J);
p=100.*[0.5/100 length(flag)/length(ecart_T) length(x)/length(I) length(y)/length(J)];

flag = find(ecart_T>0.45);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.45/100 length(flag)/length(ecart_T) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_T>0.4);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.4/100 length(flag)/length(ecart_T) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_T>0.35);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.35/100 length(flag)/length(ecart_T) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_T>0.3);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.3/100 length(flag)/length(ecart_T) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_T>0.25);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.25/100 length(flag)/length(ecart_T) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_T>0.2);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.2/100 length(flag)/length(ecart_T) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_T>0.15);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.15/100 length(flag)/length(ecart_T) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_T>0.1);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.1/100 length(flag)/length(ecart_T) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_T>0.05);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.05/100 length(flag)/length(ecart_T) length(x)/length(I) length(y)/length(J)]];

flag = find(ecart_T>0.04);
x = intersect(flag,I);
y = intersect(flag,J);
p=[p;100.*[0.04/100 length(flag)/length(ecart_T) length(x)/length(I) length(y)/length(J)]];

P_ct=array2table(p,'VariableNames',{'threshold_CT','percent_LR_tot','LR_less_10','LR_more_10'})
%% save tables 
writetable(P_sa,['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/plots/' nc_att.smru_platform_code '/comparison/' nc_att.smru_platform_code '_threshold_sa.csv']);
writetable(P_ct,['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/plots/' nc_att.smru_platform_code '/comparison/' nc_att.smru_platform_code '_threshold_ct.csv']);
writetable(P_dens,['/Users/lsiegelma/Documents/PhD/Matlab/Script_salinité/LR_HR_comparison/plots/' nc_att.smru_platform_code '/comparison/' nc_att.smru_platform_code '_threshold_dens.csv']);
%% save variables
% dlmwrite('CTLR_raw_33.txt',CTLR_brut);
% dlmwrite('SALR_raw_33.txt',SALR_brut);
% dlmwrite('DensLR_raw_33.txt',DensLR_brut);

