function compare_ctd_resolution(name_ctd,conf,plot_all_profil)
% name_ctd= identifiant de la balise
% plot_all_profil = 0 ou 1, trace tous les profils LR/HR pour comparaison
% comparaison CTD BR/CTD HR
%%
matrice=[];
matricetxt={'profil_number','latitude','longitude','rms_t','rms_s','rms_d','max diff temp','max diff sal','max diff dens'};
ctdnc=regexpdir(conf.datadir,[name_ctd '_prof.nc']);
ctd_lr=ARGO_load_qc(ctdnc{1},0);
ctd_hr=ARGO_load_qc([conf.dataprof name_ctd '_CTDHR_prof.nc'],1);
im=1;
[s, mess, messid] = mkdir(datadir_plot_hr,name_ctd);

%%

densbr=sw_dens0(ctd_lr.PSAL,ctd_lr.TEMP)-1000;
denshr=sw_dens0(ctd_hr.PSAL,ctd_hr.TEMP)-1000;

%%
list_match=[];
for ii=60:65%length(ctd_lr.JULD)
    [m,J]=min(abs(ctd_lr.JULD(ii)-ctd_hr.JULD));
    J=J;
    list_match=[list_match;J];
    profil=[];
    T_ctd=[];
    S_ctd=[];
    D_ctd=[];
    ecart_D=NaN;
    ecart_S=NaN;
    ecart_T=NaN;
    
    I=find(~isnan(ctd_lr.TEMP(:,ii)));
    if length(I)>0
        T_ctd=interp1(ctd_lr.PRES(I,ii),ctd_lr.TEMP(I,ii),ctd_hr.PRES(:,J));
        S_ctd=interp1(ctd_lr.PRES(I,ii),ctd_lr.PSAL(I,ii),ctd_hr.PRES(:,J));
        D_ctd=interp1(ctd_lr.PRES(I,ii),densbr(I,ii),ctd_hr.PRES(:,J));
        ctd_lr.TEMP_interp(:,ii)=T_ctd;
        ctd_lr.PSAL_interp(:,ii)=S_ctd;
        ctd_lr.DENS_interp(:,ii)=D_ctd;
        T=find(~isnan(D_ctd) & ~isnan(denshr(:,J)));
        ecart_D=rms(denshr(T,J)-D_ctd(T));
        if length(T)>0
            matrice(im,9)=max(abs(denshr(T,J)-D_ctd(T)));
        else
            matrice(im,9)=NaN;
        end
        T=find(~isnan(S_ctd) & ~isnan(ctd_hr.PSAL(:,J)));
        ecart_S=rms(ctd_hr.PSAL(T,J)-S_ctd(T));
        if length(T)>0
            matrice(im,8)=max(abs(ctd_hr.PSAL(T,J)-S_ctd(T)));
        else
            matrice(im,8)=NaN;
        end
        T=find(~isnan(T_ctd) & ~isnan(ctd_hr.TEMP(:,J)));
        ecart_T=rms(ctd_hr.TEMP(T,J)-T_ctd(T));
        if length(T)>0
            matrice(im,7)=max(abs(ctd_hr.TEMP(T,J)-T_ctd(T)));
        else
            matrice(im,7)=NaN;
        end
    end
    matrice(im,1)=ii;
    matrice(im,2)=ctd_lr.LATITUDE(ii);
    matrice(im,3)=ctd_lr.LONGITUDE(ii);
    matrice(im,4)=ecart_T;
    matrice(im,5)=ecart_S;
    matrice(im,6)=ecart_D;
    
    
    im=im+1;
    % calcul de la salinité corrigé par le temps de réponse
    P=ctd_hr.PRES(:,J);
    T=ctd_hr.TEMP(:,J);
    S=ctd_hr.PSAL(:,J);
    
    alpha=0.1;
    beta=0.03;
    Sc = sal_cor (alpha,beta,2,T,S,P);
    Dc = sw_dens0(Sc,T)-1000;
    
    %% figure par profil
    if plot_all_profil
        paperh=16;
        paperw=24;
        clf;
        set(gcf, ...
            'Units','characters',...
            'PaperType','A4',...
            'PaperUnits','centimeters', ...
            'PaperOrientation','portrait', ...
            'PaperPositionMode','manual', ...
            'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
        subplot(2,3,1)
        hold on
        plot(ctd_lr.TEMP(:,ii),-ctd_lr.PRES(:,ii),'b');
        plot(ctd_hr.TEMP(:,J),-ctd_hr.PRES(:,J),'g')
        xlabel('Temperature')
        ylabel('depth')
        subplot(2,3,2)
        hold on
        plot(ctd_lr.PSAL(:,ii),-ctd_lr.PRES(:,ii),'b');
        plot(ctd_hr.PSAL(:,J),-ctd_hr.PRES(:,J),'g')
        
        plot(Sc,-ctd_hr.PRES(:,J),'r')
        
        xlabel('Salinity')
        ylabel('depth')
        %      leg{1}='profil BR';
        %      leg{2}='profil HR';
        %      legend1 = legend(gca,leg);
        subplot(2,3,3)
        hold on
        plot(densbr(:,ii),-ctd_lr.PRES(:,ii),'b');
        plot(denshr(:,J),-ctd_hr.PRES(:,J),'g')
        plot(Dc,-ctd_hr.PRES(:,J),'r')
        xlabel('Density')
        ylabel('depth')
        legend({'profil BR','profil HR','profil HRC'},'Location','northeast','FontSize',6,'Position',[0.85,0.8,0.08,0.1]);
        subplot(2,3,4)
        hold on
        plot(ctd_hr.TEMP(:,J)-T_ctd,-ctd_hr.PRES(:,J));
        anno=['rms : ',num2str(ecart_T)];
        title(anno);
        %txt=annotation('textbox',[0.15 0.1 0.25 0.2],'FitHeightToText','on','String',anno );
        
        xlabel('Temperature HR - Temperature BR')
        ylabel('depth')
        subplot(2,3,5)
        hold on
        plot(ctd_hr.PSAL(:,J)-S_ctd,-ctd_hr.PRES(:,J));
        anno=['rms : ',num2str(ecart_S)];
        title(anno);
        %txt=annotation('textbox',[0.42 0.1 0.2 0.2],'FitHeightToText','on','String',anno );
        
        xlabel('Salinity HR - Salinity BR')
        ylabel('depth')
        subplot(2,3,6)
        hold on
        plot(denshr(:,J)-D_ctd,-ctd_hr.PRES(:,J));
        anno=['rms : ',num2str(ecart_D)];
        title(anno);
        %txt=annotation('textbox',[0.70 0.1 0.20 0.2],'FitHeightToText','on','String',anno );
        
        xlabel('Density HR - Density BR')
        ylabel('depth')
        
        suptitle(['Profil ',num2str(ii),' nombre de points br : ',num2str(length(I)),' seal : ',name_ctd,' beta=',num2str(beta)])
        %
        %             set(legend1,...
        %                   'Position',[0.85 0.2 0.1 0.1]);
        print('-dpng','-r500',[conf.datadir_plot_hr name_ctd '/Profil_',num2str(ii),'_',name_ctd,'_beta-',strrep(num2str(beta),'.','-')]);
    end
end


%%
I=find(isnan(matrice(:,4)));
matrice(I,:)=[];

%% Figure de comparaison
[Xi,Yi]=meshgrid(ctd_hr.JULD(list_match)-ctd_hr.JULD(list_match(1)),ctd_hr.PRES(:,1));
% compare TEMP LR/ TEMP HR
diff_temp=ctd_hr.TEMP(:,list_match)-ctd_lr.TEMP_interp;

pal=pal_bluered(100);
paperh=10;
paperw=17;
fig(2,'width',paperw,'height',paperh,'fontsize',7),clf
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
title(['individu : ' name_ctd ' Difference Temperature HR - LR']);
print('-dpng','-r500',[conf.datadir_plot_hr name_ctd '/difference_temperature_hr-lr_' name_ctd]);

% compare PSAL LR/ TEMP HR
diff_sal=ctd_hr.PSAL(:,list_match)-ctd_lr.PSAL_interp;

pal=pal_bluered(100);
paperh=10;
paperw=17;
fig(2,'width',paperw,'height',paperh,'fontsize',7),clf
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
title(['individu : ' name_ctd ' Difference Salinity HR - LR']);
print('-dpng','-r500',[conf.datadir_plot_hr name_ctd '/difference_salinity_hr-lr_' name_ctd]);

% compare PSAL LR/ TEMP HR
diff_dens=denshr(:,list_match)-ctd_lr.DENS_interp;

pal=pal_bluered(100);
paperh=10;
paperw=17;
fig(2,'width',paperw,'height',paperh,'fontsize',7),clf
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
title(['individu : ' name_ctd ' Difference Density HR - LR']);
print('-dpng','-r500',[conf.datadir_plot_hr name_ctd '/difference_density_hr-lr_' name_ctd]);
%% Tableau par zone
% Zone_tab= [];
% Zone_tabtxt={'zone','Nprof','rms_t','rms_s','rms_d'};
% lonmask=ncread([conf.maindir,'..\DATA\AVISO\mask.nc'],'NbLongitudes');
% latmask=ncread([conf.maindir,'..\DATA\AVISO\mask.nc'],'NbLatitudes');
% mask=ncread([conf.maindir,'..\DATA\AVISO\mask.nc'],'mask');
% mask_profil=NaN*(1:length(matrice))';
% mask_profil=interp2(lonmask,latmask,mask,matrice(:,3),matrice(:,2),'nearest');
% for ll=1:6
%     Zone_tab(ll,1)=ll;
%     I=find(mask_profil==ll);
%     Zone_tab(ll,2)=length(I);
%     Zone_tab(ll,3)=nanmean(matrice(I,4));
%     Zone_tab(ll,4)=nanmean(matrice(I,5));
%     Zone_tab(ll,5)=nanmean(matrice(I,6));
%
%
% end
%% histogramme des max de difference et en des rms
paperh=16;
paperw=24;
fig(1,'width',paperw,'height',paperh,'fontsize',8),clf
set(gcf, ...
    'PaperType','A4',...
    'PaperUnits','centimeters', ...
    'PaperOrientation','portrait', ...
    'PaperPositionMode','manual', ...
    'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
hist(matrice(:,7),[0:0.05:1.6])
title ( 'histogramme des differences max de temperature')
print('-dpng','-r500',[conf.datadir_plot_hr name_ctd '/hist_maxdiff_temp']);

paperh=16;
paperw=24;
fig(1,'width',paperw,'height',paperh,'fontsize',8),clf
set(gcf, ...
    'PaperType','A4',...
    'PaperUnits','centimeters', ...
    'PaperOrientation','portrait', ...
    'PaperPositionMode','manual', ...
    'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
hist(matrice(:,8),[0:0.01:0.6])
title ( 'histogramme des differences max de salinite')
print('-dpng','-r500',[conf.datadir_plot_hr name_ctd '/hist_maxdiff_sal']);

paperh=16;
paperw=24;
fig(1,'width',paperw,'height',paperh,'fontsize',8),clf
set(gcf, ...
    'PaperType','A4',...
    'PaperUnits','centimeters', ...
    'PaperOrientation','portrait', ...
    'PaperPositionMode','manual', ...
    'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
hist(matrice(:,9),[0:0.01:0.6])
title ( 'histogramme des differences max de densite')
print('-dpng','-r500',[conf.datadir_plot_hr name_ctd '/hist_maxdiff_dens']);

paperh=16;
paperw=24;
fig(1,'width',paperw,'height',paperh,'fontsize',8),clf
set(gcf, ...
    'PaperType','A4',...
    'PaperUnits','centimeters', ...
    'PaperOrientation','portrait', ...
    'PaperPositionMode','manual', ...
    'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
hist(matrice(:,4),[0:0.01:0.6])
title ( 'histogramme des rms de temperature')
print('-dpng','-r500',[conf.plotdir 'resolution_comparaison/' name_ctd '/hist_rms_temp']);

paperh=16;
paperw=24;
fig(1,'width',paperw,'height',paperh,'fontsize',8),clf
set(gcf, ...
    'PaperType','A4',...
    'PaperUnits','centimeters', ...
    'PaperOrientation','portrait', ...
    'PaperPositionMode','manual', ...
    'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
hist(matrice(:,5),[0:0.005:0.4])
title ( 'histogramme des rms de salinite')
print('-dpng','-r500',[conf.datadir_plot_hr name_ctd '/hist_rms_sal']);

paperh=16;
paperw=24;
fig(1,'width',paperw,'height',paperh,'fontsize',8),clf
set(gcf, ...
    'PaperType','A4',...
    'PaperUnits','centimeters', ...
    'PaperOrientation','portrait', ...
    'PaperPositionMode','manual', ...
    'PaperPosition',[1.5 29.7-1.5-paperh paperw paperh]);
hist(matrice(:,6),[0:0.005:0.3])
title ( 'histogramme des rms de densite')
print('-dpng','-r500',[conf.datadir_plot_hr name_ctd '/hist_rms_dens']);

%% sauvegarde de la matrice de comparaison

file=[conf.datadir_plot_hr name_ctd '/matrice.csv'];
fid=fopen(file,'w');
fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;%s;%s\n',matricetxt{1,1:9});
for ii=1:length(matrice)
    fprintf(fid,'%d;%2.5f;%2.5f;%2.7f;%2.7f;%2.7f;%2.7f;%2.7f;%2.7f\n',matrice(ii,1:9));
end
fclose(fid);

