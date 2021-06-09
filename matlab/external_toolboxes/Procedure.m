%%%%%%%%%%%%%%Correct salinity and Temperature using sal_cor.m and SPfromD0.m - May 2016
clc
clf
% clear all;
% close all;
format loose
format compact
dirplot='C:\Users\stationcalcul\GoogleDrive\MEOP-CTD\final_dataset_hr\plot\correction_sal\';
% nc_nameHR = 'C:\Users\stationcalcul\GoogleDrive\MEOP-CTD\final_dataset_hr\prof/ct96-09-13_CTDHR_prof.nc';
% nc_nameLR = 'C:\Users\stationcalcul\GoogleDrive\MEOP-CTD\final_dataset_ncARGO\IMOS\ct96\ct96-09-13_prof.nc';
MqHR=ARGO_load_qc(nc_nameHR,1);
MqLR=ARGO_load_qc(nc_nameLR,1);
nc_att=ncloadatt_struct(nc_nameHR);

%Load Elephant seal data High Resolution
lat  = ncread(nc_nameHR,'LATITUDE');
lon  = ncread(nc_nameHR,'LONGITUDE');
Temp  = ncread(nc_nameHR,'TEMP_ADJUSTED');
Sal  = ncread(nc_nameHR,'PSAL_ADJUSTED');
Den0 = sw_dens0(Sal,Temp);
Pres  = ncread(nc_nameHR,'PRES_ADJUSTED');
Date  = ncread(nc_nameHR,'JULD_LOCATION');

%Load Elephant seal data Low Resolution
latL  = ncread(nc_nameLR,'LATITUDE');
lonL  = ncread(nc_nameLR,'LONGITUDE');
TempL  = ncread(nc_nameLR,'TEMP_ADJUSTED');
SalL  = ncread(nc_nameLR,'PSAL_ADJUSTED');
Den0L = sw_dens0(SalL,TempL);
PresL  = ncread(nc_nameLR,'PRES_ADJUSTED');
DateL  = ncread(nc_nameLR,'JULD_LOCATION');

%% Procedure
mkdir([dirplot nc_att.smru_platform_code '\profils\']);
TempHR_cor=[];
SalHR_cor=[];
DensHR_cor=[];
TempHR_brut=[];
SalHR_brut=[];
DensHR_brut=[];
TempLR_brut=[];
TempLR_cor=[];
SalLR_cor=[];
SalLR_brut=[];
DensLR_cor=[];
DensLR_brut=[];
nbpt_brut=[];
d_date=[];
%plot tous les profils 1 ou 0
all_plot=0;
ind=1;
clear ind_lr d_date
for k = 1:size(Temp,2)
    [m,DP]=min(abs((Date(k)-DateL)));
    if m<0.0069 % correspond à 10min en date matlab
        %% correction pour haute resolution
        T = Temp(:,k);
        S = Sal(:,k);
        D = Den0(:,k);
        p = Pres(:,k);
        d_date(ind)=DateL(DP);
        ind_lr(ind)=DP;
        
        
        % Step 1) Interpolate S and T and get dTc using the time lag correction procedure
        %(alpha1 = .06,beta1 = 10).
        % Remove data points when dTc exceeds a certain threshold
        % (abs(dTc)<0.03), reinterpolate the density on our levels p and
        % reconstruct the salinity with this density Dc2.
        % Dc2 and Sc2 are in red on the figures.
        I = find(~isnan(p));
        Pi=(0:1000)';
        Ti=interp1(p(I),T(I),Pi);
        Si=interp1(p(I),S(I),Pi);
        [aux dTc_i] = sal_cor(.06,10,1,Ti,Si,Pi);
        dTc2 = p*NaN;
        dTc(I)=interp1(Pi,dTc_i,p(I));
        dTc2=dTc;
        
        %J=find(abs(dTc(I))<.03);
        J=1:length(I);
        Dc2 = p*NaN;
        Dc2(I) = interp1(p(I(J)),D(I(J)),p(I));
        Sc2 = SPfromD0(Dc2,T);
        
        %Step 2) Apply again the time lag correction procedure on Sc2 and the initial T (alpha2 = 0.06, beta2 = 0.04)
        % to obtain the final salinity (Sc3), the final temperature (Tc3) by adding 2*dTc, and the final density (Dc3) from Sc3 and Tcr.
        % Sc3,Tc3 and Dc3 are in blue on the figures.
        Ti=interp1(p(I),T(I),Pi);
        Si=interp1(p(I),Sc2(I),Pi);
        [Sc3_i dTc_i] = sal_cor(.06,.04,1,Ti,Si,Pi);
        dTc = p*NaN;
        dTc(I)=interp1(Pi,dTc_i,p(I));
        Tc3=T+.2*dTc;
        Sc3 = p*NaN;
        Sc3(I)=interp1(Pi,Sc3_i,p(I));
        Dc3 = sw_dens0(Sc3,Tc3);
        TempHR_cor(:,ind)=Tc3;
        SalHR_cor(:,ind)=Sc3;
        DensHR_cor(:,ind)=Dc3;
        TempHR_brut(:,ind)=T;
        SalHR_brut(:,ind)=S;
        DensHR_brut(:,ind)=D;
        %% correction pour la basse resolution
        TL = TempL(:,DP);
        SL = SalL(:,DP);
        DL = Den0L(:,DP);
        pL = PresL(:,DP);
        
        % Step 1) Interpolate S and T and get dTc using the time lag correction procedure
        %(alpha1 = .06,beta1 = 10).
        % Remove data points when dTc exceeds a certain threshold
        % (abs(dTc)<0.03), reinterpolate the density on our levels p and
        % reconstruct the salinity with this density Dc2.
        % Dc2 and Sc2 are in red on the figures.
        I = find(~isnan(pL));
        nbpt_brut(ind,1)=length(I);
        Pi=(0:1000)';
        TiL=interp1(pL(I),TL(I),Pi);
        SiL=interp1(pL(I),SL(I),Pi);
        DensLR_brut(:,ind)=interp1(pL(I),DL(I),Pi);
        TempLR_brut(:,ind)=TiL;
        SalLR_brut(:,ind)=SiL;
        [aux dTcL_i] = sal_cor(.06,10,1,TiL,SiL,Pi);
        dTc2L = pL*NaN;
        dTcL(I)=interp1(Pi,dTcL_i,pL(I));
        dTc2L=dTcL;
        
        %J=find(abs(dTc(I))<.03);
        J=1:length(I);
        Dc2L = pL*NaN;
        Dc2L(I) = interp1(pL(I(J)),DL(I(J)),pL(I));
        Sc2L = SPfromD0(Dc2L,TL);
        
        %Step 2) Apply again the time lag correction procedure on Sc2 and the initial T (alpha2 = 0.06, beta2 = 0.04)
        % to obtain the final salinity (Sc3), the final temperature (Tc3) by adding 2*dTc, and the final density (Dc3) from Sc3 and Tcr.
        % Sc3,Tc3 and Dc3 are in blue on the figures.
        TiL=interp1(pL(I),TL(I),Pi);
        SiL=interp1(pL(I),Sc2L(I),Pi);
        [Sc3L_i dTcL_i] = sal_cor(.06,.04,1,TiL,SiL,Pi);
        dTcL = pL*NaN;
        dTcL(I)=interp1(Pi,dTcL_i,pL(I));
        Tc3L=TL+.2*dTcL;
        Sc3L = pL*NaN;
        Sc3L(I)=interp1(Pi,Sc3L_i,pL(I));
        Dc3L = sw_dens0(Sc3L,Tc3L);
        TempLR_cor(:,ind)=interp1(pL(I),Tc3L(I),Pi);
        SalLR_cor(:,ind)=interp1(pL(I),Sc3L(I),Pi);
        DensLR_cor(:,ind)=interp1(pL(I),Dc3L(I),Pi);
        %% Plotting
        % plot profil HR et LR avec les 2 phases de correction
        if (all_plot & mod(k,5)==0 & sum(~isnan(T))>0)
            I=find(~isnan(T));
            subplot(2,4,1)
            plot(T,-p,'k-',Tc3,-p,'b-')
            xlabel('HR Temperature')
            xlim([min(T)-0.5 max(T)+0.5])
            subplot(2,4,2)
            plot(S,-p,'k-',Sc2,-p,'r-',Sc3,-p,'b-')
            xlabel('HR Salinity')
            xlim([min(S)-0.3 max(S)+0.3])
            subplot(2,4,3)
            plot(D,-p,'k-',Dc2,-p,'r-',Dc3,-p,'b-')
            xlabel('HR Sigma 0')
            xlim([1026 1028])
            %         legend('original profile','profile step 1','profile step 2','Location','southeast')
            subplot(2,4,4)
            plot(dTc(I),-p(I),'k-',dTc2(I),-p(I),'r-')
            xlabel('HR dTc')
            xlim([min(dTc)-0.1 max(dTc)+0.1])
            I=find(~isnan(TL));
            subplot(2,4,5)
            plot(TL,-pL,'k-',Tc3L,-pL,'b-')
            xlabel('LR Temperature')
            xlim([min(T)-0.5 max(T)+0.5])
            subplot(2,4,6)
            plot(SL,-pL,'k-',Sc2L,-pL,'r-',Sc3L,-pL,'b-')
            xlabel('LR Salinity')
            xlim([min(S)-0.3 max(S)+0.3])
            subplot(2,4,7)
            plot(DL,-pL,'k-',Dc2L,-pL,'r-',Dc3L,-pL,'b-')
            xlabel('LR Sigma 0')
            xlim([1026 1028])
            subplot(2,4,8)
            plot(dTcL(I),-pL(I),'k-',dTc2L(I),-pL(I),'r-')
            xlabel('LR dTc')
            xlim([min(dTc)-0.1 max(dTc)+0.1])
            print('-dpng',[dirplot nc_att.smru_platform_code '\profils\' nc_att.smru_platform_code '_profilHR_' num2str(k) '.png'],'-r0');
        end
        ind=ind+1;
    end
end

%% plot comparaison des corrections HR et LR
if all_plot
    ecart_T=[];
    ecart_S=[];
    ecart_D=[];
    mkdir([dirplot nc_att.smru_platform_code '\comparaison\profils\']);
    for ii=1:10:size(TempHR_cor,2)
        subplot(2,3,1)
        plot(TempHR_brut(:,ii),-p,'k-',TempHR_cor(:,ii),-p,'b-',TempLR_cor(:,ii),-p,'r-')
        xlabel('Temperature')
        if sum(~isnan(TempHR_cor(:,ii)))
            xlim([min(TempHR_cor(:,ii))-0.5 max(TempHR_cor(:,ii))+0.5])
        else
            xlim([-0.5 5])
        end
        subplot(2,3,2)
        plot(SalHR_brut(:,ii),-p,'k-',SalHR_cor(:,ii),-p,'b-',SalLR_cor(:,ii),-p,'r-')
        xlabel('Salinity')
        if sum(~isnan(SalHR_cor(:,ii)))
            xlim([min(SalHR_cor(:,ii))-0.3 max(SalHR_cor(:,ii))+0.3])
        else
            xlim([33 35])
        end
        
        subplot(2,3,3)
        plot(DensHR_brut(:,ii),-p,'k-',DensHR_cor(:,ii),-p,'b-',DensLR_cor(:,ii),-p,'r-')
        xlabel('HR Sigma 0')
        xlim([1026 1028])
        legend('HR raw','HR corrected','LR corrected','Location','southeast')
        
        subplot(2,3,4)
        plot(TempHR_cor(:,ii)-TempLR_cor(:,ii),-p,'k-')
        T=find(~isnan(TempLR_cor(:,ii)) & ~isnan(TempHR_cor(:,ii)));
        ecart_T(ii,1)=rms(TempHR_cor(T,ii)-TempLR_cor(T,ii));
        anno=['rms : ',num2str(ecart_T(ii,1))];
        title(anno);
        xlabel('THR-TLR')
        %xlim([min(T)-0.5 max(T)+0.5])
        
        subplot(2,3,5)
        plot(SalHR_cor(:,ii)-SalLR_cor(:,ii),-p,'k-')
        T=find(~isnan(SalLR_cor(:,ii)) & ~isnan(SalHR_cor(:,ii)));
        ecart_S(ii,1)=rms(SalHR_cor(T,ii)-SalLR_cor(T,ii));
        anno=['rms : ',num2str(ecart_S(ii,1))];
        title(anno);
        xlabel('SHR-SLR')
        
        subplot(2,3,6)
        plot(DensHR_cor(:,ii)-DensLR_cor(:,ii),-p,'k-')
        T=find(~isnan(DensLR_cor(:,ii)) & ~isnan(DensHR_cor(:,ii)));
        ecart_D(ii,1)=rms(DensHR_cor(T,ii)-DensLR_cor(T,ii));
        anno=['rms : ',num2str(ecart_D(ii,1))];
        title(anno);
        xlabel('DHR-DLR')
        
        print('-dpng',[dirplot nc_att.smru_platform_code '\comparaison\profils\' nc_att.smru_platform_code '_profil_' num2str(ii) '.png'],'-r0');
        
    end
end
% fid=fopen([dirplot nc_att.smru_platform_code '\comparaison\' nc_att.smru_platform_code 'mean_rms.txt'],'w');
% fprintf(fid,'mean rms temperature : %f \n',nanmean(ecart_T));
% fprintf(fid,'mean rms salinity : %f \n',nanmean(ecart_S));
% fprintf(fid,'mean rms sigma0 : %f \n',nanmean(ecart_D));
%
% fclose(fid);



%% comparaison HR cor/ LR cor
if all_plot
    clf
    [Xi,Yi]=meshgrid(d_date-d_date(1),p);
    moy=[];
    ecart=[];
    % compare TEMP LR/ TEMP HR
    diff_temp=TempHR_cor-TempLR_cor;
    moy(:,1)=nanmean(diff_temp');
    ecart(:,1)=nanstd(diff_temp');
    pal=pal_bluered(100);
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,-Yi,diff_temp)
    shading('interp')
    colormap(pal)
    caxis([-0.4 0.4])
    colorbar
    axis([min(Xi(:))-5 max(Xi(:)) -1000 0]);
    xlabel('Days')
    ylabel('Depth')
    title(['individu : ' nc_att.smru_platform_code ' Difference Temperature Corrected HRcor - LRcor']);
    print('-dpng','-r300',[dirplot nc_att.smru_platform_code '\comparaison\' '/difference_temperature_hrcor-lrcor_' nc_att.smru_platform_code]);
    
    % compare SAL LR/ SAL HR
    diff_sal=SalHR_cor-SalLR_cor;
    moy(:,2)=nanmean(diff_sal');
    ecart(:,2)=nanstd(diff_sal');
    
    pal=pal_bluered(100);
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,-Yi,diff_sal)
    shading('interp')
    colormap(pal)
    caxis([-0.2 0.2])
    colorbar
    axis([min(Xi(:))-5 max(Xi(:)) -1000 0]);
    xlabel('Days')
    ylabel('Depth')
    title(['individu : ' nc_att.smru_platform_code ' Difference Salinity Corrected HRcor - LRcor']);
    print('-dpng','-r300',[dirplot nc_att.smru_platform_code '\comparaison\' '/difference_salinity_hrcor-lrcor_' nc_att.smru_platform_code]);
    
    % compare DENS LR/ DENS HR
    diff_dens=DensHR_cor-DensLR_cor;
    moy(:,3)=nanmean(diff_dens');
    ecart(:,3)=nanstd(diff_dens');
    
    pal=pal_bluered(100);
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,-Yi,diff_dens)
    shading('interp')
    colormap(pal)
    caxis([-0.1 0.1])
    colorbar
    axis([min(Xi(:))-5 max(Xi(:)) -1000 0]);
    xlabel('Days')
    ylabel('Depth')
    title(['individu : ' nc_att.smru_platform_code ' Difference Density Corrected HRcor - LRcor']);
    print('-dpng','-r300',[dirplot nc_att.smru_platform_code '\comparaison\' '/difference_density_hrcor-lrcor_' nc_att.smru_platform_code]);
end

%% difference entre HR raw et HR corrected
if all_plot
    clf
    [Xi,Yi]=meshgrid(d_date-d_date(1),p);
    moy=[];
    ecart=[];
    % compare TEMP HR/ TEMP HR cor
    diff_temp=TempHR_cor-TempHR_brut;
    moy(:,1)=nanmean(diff_temp');
    ecart(:,1)=nanstd(diff_temp');
    pal=pal_bluered(100);
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,-Yi,diff_temp)
    shading('interp')
    colormap(pal)
    caxis([-0.4 0.4])
    colorbar
    axis([min(Xi(:))-5 max(Xi(:)) -1000 0]);
    xlabel('Days')
    ylabel('Depth')
    title(['individu : ' nc_att.smru_platform_code ' Difference Temperature Corrected HR - Raw HR']);
    print('-dpng','-r300',[dirplot nc_att.smru_platform_code '\comparaison\' '/difference_temperature_hrcor-hrraw_' nc_att.smru_platform_code]);
    
    % compare SAL HR/ SAL HR raw
    diff_sal=SalHR_cor-SalHR_brut;
    moy(:,2)=nanmean(diff_sal');
    ecart(:,2)=nanstd(diff_sal');
    
    pal=pal_bluered(100);
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,-Yi,diff_sal)
    shading('interp')
    colormap(pal)
    caxis([-0.2 0.2])
    colorbar
    axis([min(Xi(:))-5 max(Xi(:)) -1000 0]);
    xlabel('Days')
    ylabel('Depth')
    title(['individu : ' nc_att.smru_platform_code ' Difference Salinity Corrected HR - HR raw']);
    print('-dpng','-r300',[dirplot nc_att.smru_platform_code '\comparaison\' '/difference_salinity_hrcor-hrraw_' nc_att.smru_platform_code]);
    
    % compare DENS HR/ DENS HR raw
    diff_dens=DensHR_cor-DensHR_brut;
    moy(:,3)=nanmean(diff_dens');
    ecart(:,3)=nanstd(diff_dens');
    
    pal=pal_bluered(100);
    paperh=10;
    paperw=17;
    set(gcf, ...
        'PaperType','A4',...
        'PaperUnits','centimeters', ...
        'PaperOrientation','portrait', ...
        'PaperPositionMode','manual', ...
        'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
    hold on
    pcolor(Xi,-Yi,diff_dens)
    shading('interp')
    colormap(pal)
    caxis([-0.1 0.1])
    colorbar
    axis([min(Xi(:))-5 max(Xi(:)) -1000 0]);
    xlabel('Days')
    ylabel('Depth')
    title(['individu : ' nc_att.smru_platform_code ' Difference Density Corrected HR - HR raw']);
    print('-dpng','-r300',[dirplot nc_att.smru_platform_code '\comparaison\' '/difference_density_hrcor-hrraw_' nc_att.smru_platform_code]);
    
end
%% comparaison des profils basse resolution en fonction du nombre de points par profils
clf
% on classe les profils en fonction du nombre de points
% seuil à 10 points
I=find(sum(~isnan(MqLR.PSAL(:,ind_lr)))<10);
J=find(sum(~isnan(MqLR.PSAL(:,ind_lr)))>=10);
diff_temp=TempHR_cor-TempLR_cor;
diff_sal=SalHR_cor-SalLR_cor;
diff_dens=DensHR_cor-DensLR_cor;
I_all=[I_all;I'+length(temp_HR_LR)];
J_all=[J_all;J'+length(temp_HR_LR)];

temp_HR_LR=[temp_HR_LR,diff_temp];
sal_HR_LR=[sal_HR_LR,diff_sal];
dens_HR_LR=[dens_HR_LR,diff_dens];
% plot des differences mean et std HR cor et LR cor
% LR divisé en 3 :all ,% LR-10pts , LR+10pts


subplot(2,3,1)
plot(nanmean(diff_temp'),-p,'-k',nanmean(diff_temp(:,I)'),-p,'-b',nanmean(diff_temp(:,J)'),-p,'-r')
xlabel('Mean diff temp')
ylabel('Depth')
subplot(2,3,2)
plot(nanmean(diff_sal'),-p,'-k',nanmean(diff_sal(:,I)'),-p,'-b',nanmean(diff_sal(:,J)'),-p,'-r')
xlabel('Mean diff sal')
subplot(2,3,3)
plot(nanmean(diff_dens'),-p,'-k',nanmean(diff_dens(:,I)'),-p,'-b',nanmean(diff_dens(:,J)'),-p,'-r')
legend('all','LR -10pts','LR +10pts','Location','southwest','FontSize',2)
xlabel('Mean diff sigma0')
subplot(2,3,4)
plot(nanstd(diff_temp'),-p,'-k',nanstd(diff_temp(:,I)'),-p,'-b',nanstd(diff_temp(:,J)'),-p,'-r')
xlabel('STD diff temp')
ylabel('Depth')
subplot(2,3,5)
plot(nanstd(diff_sal'),-p,'-k',nanstd(diff_sal(:,I)'),-p,'-b',nanstd(diff_sal(:,J)'),-p,'-r')
xlabel('STD diff temp')
subplot(2,3,6)
plot(nanstd(diff_dens'),-p,'-k',nanstd(diff_dens(:,I)'),-p,'-b',nanstd(diff_dens(:,J)'),-p,'-r')
xlabel('STD diff temp')
suptitle(['individu : ' nc_att.smru_platform_code ' Profil mean and std difference HR-LR with number points LR threshold = 10']);
print([dirplot nc_att.smru_platform_code '\comparaison\' '/mean_std_diff_nbpts_LR_' nc_att.smru_platform_code],'-depsc');
print([dirplot nc_att.smru_platform_code '\comparaison\' '/mean_std_diff_nbpts_LR_' nc_att.smru_platform_code],'-dpng','-r400');

%% mean/std difference LAG et raw
clf
I=find(sum(~isnan(MqLR.PSAL(:,ind_lr)))<10);
J=find(sum(~isnan(MqLR.PSAL(:,ind_lr)))>=10);
diff_temp=TempLR_cor-TempLR_brut;
diff_sal=SalLR_cor-SalLR_brut;
diff_dens=DensLR_cor-DensLR_brut;
diff_tempHR=TempHR_cor-TempHR_brut;
diff_salHR=SalHR_cor-SalHR_brut;
diff_densHR=DensHR_cor-DensHR_brut;
diff_temp_LR=[diff_temp_LR,diff_temp];
diff_sal_LR=[diff_sal_LR,diff_sal];
diff_dens_LR=[diff_dens_LR,diff_dens];
diff_temp_HR=[diff_temp_HR,diff_tempHR];
diff_sal_HR=[diff_sal_HR,diff_salHR];
diff_dens_HR=[diff_dens_HR,diff_densHR];
% plot des differences mean et std entre la LR-10pts et HR et la LR+10pts
% et HR

subplot(2,3,1)
plot(nanmean(diff_tempHR'),-p,'-g',nanmean(diff_temp'),-p,'-k',nanmean(diff_temp(:,I)'),-p,'-b',nanmean(diff_temp(:,J)'),-p,'-r')
xlabel('Mean diff temp')
ylabel('Depth')
subplot(2,3,2)
plot(nanmean(diff_salHR'),-p,'-g',nanmean(diff_sal'),-p,'-k',nanmean(diff_sal(:,I)'),-p,'-b',nanmean(diff_sal(:,J)'),-p,'-r')
xlabel('Mean diff sal')
subplot(2,3,3)
plot(nanmean(diff_densHR'),-p,'-g',nanmean(diff_dens'),-p,'-k',nanmean(diff_dens(:,I)'),-p,'-b',nanmean(diff_dens(:,J)'),-p,'-r')
legend('HR','all LR','LR -10pts','LR +10pts','Location','southwest','FontSize',2)
xlabel('Mean diff sigma0')
subplot(2,3,4)
plot(nanstd(diff_tempHR'),-p,'-g',nanstd(diff_temp'),-p,'-k',nanstd(diff_temp(:,I)'),-p,'-b',nanstd(diff_temp(:,J)'),-p,'-r')
xlabel('STD diff temp')
ylabel('Depth')
subplot(2,3,5)
plot(nanstd(diff_salHR'),-p,'-g',nanstd(diff_sal'),-p,'-k',nanstd(diff_sal(:,I)'),-p,'-b',nanstd(diff_sal(:,J)'),-p,'-r')
xlabel('STD diff temp')
subplot(2,3,6)
plot(nanstd(diff_densHR'),-p,'-g',nanstd(diff_dens'),-p,'-k',nanstd(diff_dens(:,I)'),-p,'-b',nanstd(diff_dens(:,J)'),-p,'-r')
xlabel('STD diff temp')
suptitle(['individu : ' nc_att.smru_platform_code ' Profil mean and std difference cor - raw']);
print([dirplot nc_att.smru_platform_code '\comparaison\' '/mean_std_diff_cor-raw_' nc_att.smru_platform_code],'-depsc');
print([dirplot nc_att.smru_platform_code '\comparaison\' '/mean_std_diff_cor-raw_' nc_att.smru_platform_code],'-dpng','-r400');



%% smooth des profils HR
clf
mkdir([dirplot nc_att.smru_platform_code '\comparaison\smooth\']);
mkdir([dirplot nc_att.smru_platform_code '\comparaison\raw_cal\']);
mkdir([dirplot nc_att.smru_platform_code '\comparaison\resolution\']);
TempHR_sm=NaN*TempHR_cor;
SalHR_sm=NaN*SalHR_cor;
DensHR_sm=NaN*DensHR_cor;

for ii=1:size(TempHR_cor,2)
    I=find(~isnan(TempHR_cor(:,ii)));
    %TempHR_sm(I,ii)=smooth(TempHR_cor(I,ii),0.15,'loess');
    TempHR_sm(I,ii)=smooth(TempHR_cor(I,ii),11,'moving');
    I=find(~isnan(SalHR_cor(:,ii)));
    SalHR_sm(I,ii)=smooth(SalHR_cor(I,ii),11,'moving');
    
end
DensHR_sm = sw_dens0(SalHR_sm,TempHR_sm);
if all_plot
    for ii=1:10:size(TempHR_cor,2)
        %% figure HR raw/cor et smooth
        clf
        subplot(2,3,1)
        plot(TempHR_brut(:,ii),-p,'k-',TempHR_cor(:,ii),-p,'b-',TempHR_sm(:,ii),-p,'r-')
        xlabel('Temperature')
        if sum(~isnan(TempHR_cor(:,ii)))
            xlim([min(TempHR_cor(:,ii))-0.5 max(TempHR_cor(:,ii))+0.5])
        else
            xlim([-0.5 5])
        end
        
        subplot(2,3,2)
        plot(SalHR_brut(:,ii),-p,'k-',SalHR_cor(:,ii),-p,'b-',SalHR_sm(:,ii),-p,'r-')
        xlabel('Salinity')
        if sum(~isnan(SalHR_cor(:,ii)))
            xlim([min(SalHR_cor(:,ii))-0.3 max(SalHR_cor(:,ii))+0.3])
        else
            xlim([33 35])
        end
        subplot(2,3,3)
        plot(DensHR_brut(:,ii),-p,'k-',DensHR_cor(:,ii),-p,'b-',DensHR_sm(:,ii),-p,'r-')
        xlabel('HR Sigma 0')
        xlim([1026 1028])
        legend('HR raw','HR lag','HR lag+smooth','Location','southeast')
        
        subplot(2,3,4)
        plot(TempHR_cor(:,ii)-TempHR_sm(:,ii),-p,'k-')
        T=find(~isnan(TempHR_sm(:,ii)) & ~isnan(TempHR_cor(:,ii)));
        ecart_T(ii,1)=rms(TempHR_cor(T,ii)-TempHR_sm(T,ii));
        anno=['rms : ',num2str(ecart_T(ii,1))];
        title(anno);
        xlabel('THR-THRsmooth')
        %xlim([min(T)-0.5 max(T)+0.5])
        
        subplot(2,3,5)
        plot(SalHR_cor(:,ii)-SalHR_sm(:,ii),-p,'k-')
        T=find(~isnan(SalHR_sm(:,ii)) & ~isnan(SalHR_cor(:,ii)));
        ecart_S(ii,1)=rms(SalHR_cor(T,ii)-SalHR_sm(T,ii));
        anno=['rms : ',num2str(ecart_S(ii,1))];
        title(anno);
        xlabel('SHR-SHRsmooth')
        
        subplot(2,3,6)
        plot(DensHR_cor(:,ii)-DensHR_sm(:,ii),-p,'k-')
        T=find(~isnan(DensHR_sm(:,ii)) & ~isnan(DensHR_cor(:,ii)));
        ecart_D(ii,1)=rms(DensHR_cor(T,ii)-DensHR_sm(T,ii));
        anno=['rms : ',num2str(ecart_D(ii,1))];
        title(anno);
        xlabel('DHR-DHRsmooth')
        suptitle(['individu : ' nc_att.smru_platform_code ' Profil'  num2str(ii) ' HR smooth method : moving and span = 11']);
        print('-dpng',[dirplot nc_att.smru_platform_code '\comparaison\smooth\' nc_att.smru_platform_code '_profil_' num2str(ii) '_smooth_moving.png'],'-r400');
        
        %% figure qui trace les profils HR raw/Lag/smooth et LR raw/LAG
        clf
        subplot(2,3,1)
        plot(TempHR_brut(:,ii),-p,'k-',TempLR_brut(:,ii),-p,'b-',TempHR_sm(:,ii),-p,'r-',TempLR_cor(:,ii),-p,'g-')
        xlabel('Temperature')
        if sum(~isnan(TempHR_cor(:,ii)))
            xlim([min(TempHR_cor(:,ii))-0.5 max(TempHR_cor(:,ii))+0.5])
        else
            xlim([-0.5 5])
        end
        
        subplot(2,3,2)
        plot(SalHR_brut(:,ii),-p,'k-',SalLR_brut(:,ii),-p,'b-',SalHR_sm(:,ii),-p,'r-',SalLR_cor(:,ii),-p,'g-')
        xlabel('Salinity')
        if sum(~isnan(SalHR_cor(:,ii)))
            xlim([min(SalHR_cor(:,ii))-0.3 max(SalHR_cor(:,ii))+0.3])
        else
            xlim([33 35])
        end
        subplot(2,3,3)
        plot(DensHR_brut(:,ii),-p,'k-',DensLR_brut(:,ii),-p,'b-',DensHR_sm(:,ii),-p,'r-',DensLR_cor(:,ii),-p,'g-')
        xlabel('Sigma 0')
        xlim([1026 1028])
        legend('HR raw','LR raw','HR lag+smooth','LR lag','Location','southeast')
        
        subplot(2,3,4)
        plot(TempHR_brut(:,ii)-TempHR_sm(:,ii),-p,'b-',TempLR_brut(:,ii)-TempLR_cor(:,ii),-p,'r-')
        T=find(~isnan(TempHR_sm(:,ii)) & ~isnan(TempHR_cor(:,ii)));
        ecart_T(ii,1)=rms(TempHR_brut(T,ii)-TempHR_sm(T,ii));
        T=find(~isnan(TempLR_brut(:,ii)) & ~isnan(TempLR_cor(:,ii)));
        anno=['rms : ',num2str(ecart_T(ii,1),'%.3f'),' rms LR : ',num2str(rms(TempLR_brut(T,ii)-TempLR_cor(T,ii)),'%.3f')];
        title(anno);
        xlabel('T-Tcor')
        %xlim([min(T)-0.5 max(T)+0.5])
        
        subplot(2,3,5)
        plot(SalHR_brut(:,ii)-SalHR_sm(:,ii),-p,'b-',SalLR_brut(:,ii)-SalLR_cor(:,ii),-p,'r-')
        T=find(~isnan(SalHR_sm(:,ii)) & ~isnan(SalHR_cor(:,ii)));
        ecart_S(ii,1)=rms(SalHR_brut(T,ii)-SalHR_sm(T,ii));
        T=find(~isnan(SalLR_brut(:,ii)) & ~isnan(SalLR_cor(:,ii)));
        anno=['rms : ',num2str(ecart_S(ii,1),'%.3f'),' rms LR : ',num2str(rms(SalLR_brut(T,ii)-SalLR_cor(T,ii)),'%.3f')];
        title(anno);
        xlabel('S-Scor')
        
        subplot(2,3,6)
        plot(DensHR_brut(:,ii)-DensHR_sm(:,ii),-p,'b-',DensLR_brut(:,ii)-DensLR_cor(:,ii),-p,'r-')
        T=find(~isnan(DensHR_sm(:,ii)) & ~isnan(DensHR_cor(:,ii)));
        ecart_D(ii,1)=rms(DensHR_brut(T,ii)-DensHR_sm(T,ii));
        T=find(~isnan(DensLR_brut(:,ii)) & ~isnan(DensLR_cor(:,ii)));
        anno=['rms HR : ',num2str(ecart_D(ii,1),'%.3f'),' rms LR : ',num2str(rms(DensLR_brut(T,ii)-DensLR_cor(T,ii)),'%.3f')];
        legend('HR raw - HR Lagsmooth','LR raw - LR Lag','Location','southeast')
        title(anno);
        xlabel('D-Dcor')
        suptitle(['individu : ' nc_att.smru_platform_code ' Profil '  num2str(ii) ' nb points LR : ' num2str(nbpt_brut(ii))]);
        print('-dpng',[dirplot nc_att.smru_platform_code '\comparaison\raw_cal\' nc_att.smru_platform_code '_profil_' num2str(ii) '.png'],'-r400');
        
        %% figure comparaison HR Lag smooth et LR Lag
        clf
        subplot(2,3,1)
        plot(TempHR_sm(:,ii),-p,'b-',TempLR_cor(:,ii),-p,'r-')
        xlabel('Temperature')
        if sum(~isnan(TempHR_cor(:,ii)))
            xlim([min(TempHR_cor(:,ii))-0.5 max(TempHR_cor(:,ii))+0.5])
        else
            xlim([-0.5 5])
        end
        
        subplot(2,3,2)
        plot(SalHR_sm(:,ii),-p,'b-',SalLR_cor(:,ii),-p,'r-')
        xlabel('Salinity')
        if sum(~isnan(SalHR_cor(:,ii)))
            xlim([min(SalHR_cor(:,ii))-0.3 max(SalHR_cor(:,ii))+0.3])
        else
            xlim([33 35])
        end
        subplot(2,3,3)
        plot(DensHR_sm(:,ii),-p,'b-',DensLR_cor(:,ii),-p,'r-')
        xlabel('Sigma 0')
        xlim([1026 1028])
        legend('HR lag+smooth','LR lag','Location','southeast')
        
        subplot(2,3,4)
        plot(TempHR_sm(:,ii)-TempLR_cor(:,ii),-p,'b-')
        T=find(~isnan(TempHR_sm(:,ii)) & ~isnan(TempLR_cor(:,ii)));
        ecart_T(ii,1)=rms(TempHR_sm(T,ii)-TempLR_cor(T,ii));
        anno=['rms : ',num2str(ecart_T(ii,1))];
        title(anno);
        xlabel('THR-TLR')
        
        subplot(2,3,5)
        plot(SalHR_sm(:,ii)-SalLR_cor(:,ii),-p,'b-')
        T=find(~isnan(SalHR_sm(:,ii)) & ~isnan(SalLR_cor(:,ii)));
        ecart_S(ii,1)=rms(SalHR_sm(T,ii)-SalLR_cor(T,ii));
        anno=['rms : ',num2str(ecart_S(ii,1))];
        title(anno);
        xlabel('SHR-SLR')
        
        subplot(2,3,6)
        plot(DensHR_sm(:,ii)-DensLR_cor(:,ii),-p,'b-')
        T=find(~isnan(DensHR_sm(:,ii)) & ~isnan(DensLR_cor(:,ii)));
        ecart_D(ii,1)=rms(DensHR_sm(T,ii)-DensLR_cor(T,ii));
        anno=['rms HR : ',num2str(ecart_D(ii,1))];
        legend('HR Lagsmooth - LR Lag','Location','southeast')
        title(anno);
        xlabel('DHR-DLR')
        suptitle(['individu : ' nc_att.smru_platform_code ' Profil '  num2str(ii) ' nb points LR : ' num2str(nbpt_brut(ii))]);
        print('-dpng',[dirplot nc_att.smru_platform_code '\comparaison\resolution\' nc_att.smru_platform_code '_profil_' num2str(ii) '.png'],'-r400');
        
    end
end
%% figure des transec T/S/D de 0 à 500m
if all_plot
    colormap default
    Tcontour=[-1:2:6];
    Scontour=[33.6:0.5:34.6];
    Dcontour=[1026.5:0.5:1027.5];
    pmax=500;
    [Xi,Yi]=meshgrid(d_date-d_date(1),p(1:pmax));
    
    % figure temperature HR raw
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
    pcolor(Xi,-Yi,TempHR_brut(1:pmax,:))
    shading('interp')
    contour(Xi,-Yi,TempHR_brut(1:pmax,:),Tcontour,'linecolor','k')
    caxis([-1 6])
    colorbar('northoutside')
    title('Temperature HR raw');
    print('-dpng',[dirplot nc_att.smru_platform_code '\comparaison\' nc_att.smru_platform_code '_transec_tempHR_raw.png'],'-r300');
    
    % figure temperature HR cor
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
    pcolor(Xi,-Yi,TempHR_cor(1:pmax,:))
    shading('interp')
    contour(Xi,-Yi,TempHR_cor(1:pmax,:),Tcontour,'linecolor','k')
    caxis([-1 6])
    colorbar('northoutside')
    title('Temperature HR cor');
    print('-dpng',[dirplot nc_att.smru_platform_code '\comparaison\' nc_att.smru_platform_code '_transec_tempHR_cor.png'],'-r300');
    
    % figure temperature LR raw
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
    pcolor(Xi,-Yi,TempLR_brut(1:pmax,:))
    shading('interp')
    contour(Xi,-Yi,TempLR_brut(1:pmax,:),Tcontour,'linecolor','k')
    caxis([-1 6])
    colorbar('northoutside')
    title('Temperature LR raw');
    print('-dpng',[dirplot nc_att.smru_platform_code '\comparaison\' nc_att.smru_platform_code '_transec_tempLR_raw.png'],'-r300');
    
    % figure temperature LR cor
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
    pcolor(Xi,-Yi,TempLR_cor(1:pmax,:))
    shading('interp')
    contour(Xi,-Yi,TempLR_cor(1:pmax,:),Tcontour,'linecolor','k')
    caxis([-1 6])
    colorbar('northoutside')
    title('Temperature LR cor');
    print('-dpng',[dirplot nc_att.smru_platform_code '\comparaison\' nc_att.smru_platform_code '_transec_tempLR_cor.png'],'-r300');
    
    % figure salinity HR raw
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
    pcolor(Xi,-Yi,SalHR_brut(1:pmax,:))
    shading('interp')
    contour(Xi,-Yi,SalHR_brut(1:pmax,:),Scontour,'linecolor','k')
    caxis([33.6 34.6])
    colorbar('northoutside')
    title('Salinity HR raw');
    print('-dpng',[dirplot nc_att.smru_platform_code '\comparaison\' nc_att.smru_platform_code '_transec_salHR_raw.png'],'-r300');
    
    % figure salinity HR cor
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
    pcolor(Xi,-Yi,SalHR_cor(1:pmax,:))
    shading('interp')
    contour(Xi,-Yi,SalHR_cor(1:pmax,:),Scontour,'linecolor','k')
    caxis([33.6 34.6])
    colorbar('northoutside')
    title('Salinity HR cor');
    print('-dpng',[dirplot nc_att.smru_platform_code '\comparaison\' nc_att.smru_platform_code '_transec_salHR_cor.png'],'-r300');
    
    % figure salinity LR raw
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
    pcolor(Xi,-Yi,SalLR_brut(1:pmax,:))
    shading('interp')
    contour(Xi,-Yi,SalLR_brut(1:pmax,:),Scontour,'linecolor','k')
    caxis([33.6 34.6])
    colorbar('northoutside')
    title('Salinity LR raw');
    print('-dpng',[dirplot nc_att.smru_platform_code '\comparaison\' nc_att.smru_platform_code '_transec_salLR_raw.png'],'-r300');
    
    % figure salinite LR cor
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
    pcolor(Xi,-Yi,SalLR_cor(1:pmax,:))
    shading('interp')
    contour(Xi,-Yi,SalLR_cor(1:pmax,:),Scontour,'linecolor','k')
    caxis([33.6 34.6])
    colorbar('northoutside')
    title('Salinity LR cor');
    print('-dpng',[dirplot nc_att.smru_platform_code '\comparaison\' nc_att.smru_platform_code '_transec_salLR_cor.png'],'-r300');
    
    % figure densite HR raw
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
    pcolor(Xi,-Yi,DensHR_brut(1:pmax,:))
    shading('interp')
    contour(Xi,-Yi,DensHR_brut(1:pmax,:),Dcontour,'linecolor','k')
    caxis([1026.5 1027.5])
    colorbar('northoutside')
    title('Sigma0 HR raw');
    print('-dpng',[dirplot nc_att.smru_platform_code '\comparaison\' nc_att.smru_platform_code '_transec_densHR_raw.png'],'-r300');
    
    % figure densite HR cor
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
    pcolor(Xi,-Yi,DensHR_cor(1:pmax,:))
    shading('interp')
    contour(Xi,-Yi,DensHR_cor(1:pmax,:),Dcontour,'linecolor','k')
    caxis([1026.5 1027.5])
    colorbar('northoutside')
    title('Sigma0 HR cor');
    print('-dpng',[dirplot nc_att.smru_platform_code '\comparaison\' nc_att.smru_platform_code '_transec_densHR_cor.png'],'-r300');
    
    % figure densite LR raw
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
    pcolor(Xi,-Yi,DensLR_brut(1:pmax,:))
    shading('interp')
    contour(Xi,-Yi,DensLR_brut(1:pmax,:),Dcontour,'linecolor','k')
    caxis([1026.5 1027.5])
    colorbar('northoutside')
    title('Sigma0 LR raw');
    print('-dpng',[dirplot nc_att.smru_platform_code '\comparaison\' nc_att.smru_platform_code '_transec_densLR_raw.png'],'-r300');
    
    % figure densite LR cor
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
    pcolor(Xi,-Yi,DensLR_cor(1:pmax,:))
    shading('interp')
    contour(Xi,-Yi,DensLR_cor(1:pmax,:),Dcontour,'linecolor','k')
    caxis([1026.5 1027.5])
    colorbar('northoutside')
    title('Sigma0 LR cor');
    print('-dpng',[dirplot nc_att.smru_platform_code '\comparaison\' nc_att.smru_platform_code '_transec_densLR_cor.png'],'-r300');
end