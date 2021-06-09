function transect_time_ARGO_CHLA(argo_qc,Tcontour,Fcontour,Ccontour,hfig,printname,plot_conf)
% function transect_time_ARGO
%   plot time-referenced transects of T and CHLA
%
% INPUT:
% argo_qc: struct variable obtained using the function ARGO_load_qc
% Tcontour, Scontour: list of contours, or number of contours to use.
% Ccontour: list of contours for the date colorbar. at least two-element
%       vector. used to scale X-axis.
% hfig: figure handle number. if 0, a new figure is opened.
% printname: print name. no print if empty string. suffix indicates the print
%       format.

if argo_qc.np==0, return, end
warning off

N=argo_qc.np;
zp=argo_qc.PRES';
zt=argo_qc.TEMP';
zs=argo_qc.CHLA';

prof_lat=argo_qc.LATITUDE;
prof_lon=argo_qc.LONGITUDE;
time=argo_qc.JULD_LOCATION-min(argo_qc.JULD_LOCATION);
Pbot=zeros(N,1); 
for kk=1:N, 
    I=find(~isnan(zt(kk,:))); 
    if length(I), Pbot(kk)=zp(kk,I(end)); end
end

%% transect: data preparation

corr=1; % smoothing coefficient (unit: day)
pmax=200;
if exist('plot_conf') & isfield(plot_conf,'transect_corr'), corr=plot_conf.transect_corr; end
if exist('plot_conf') & isfield(plot_conf,'pmax_fluo'), pmax=plot_conf.pmax_fluo; end

pres = (5:5:pmax)';
zti = zeros(length(pres),N)*NaN;
zsi = zeros(length(pres),N)*NaN;
for kk=1:N,
    try
    zti(:,kk) = interp1(zp(kk,~isnan(zt(kk,:))),zt(kk,~isnan(zt(kk,:))),pres);
    zsi(:,kk) = interp1(zp(kk,~isnan(zs(kk,:))),zs(kk,~isnan(zs(kk,:))),pres);
    end
end

% lissage temp�rature
Ld=time;Pd=pres;Xd=zti;
Ld2=(1:corr:floor(max(Ld)))';Pd2=pres;Xd2=NaN*zeros(length(Pd2),length(Ld2));
for nt=1:length(Ld2),
    K=find(abs(Ld-Ld2(nt))<5*corr);
    if isempty(K), continue, end
    Ldi=ones(size(Pd))*Ld(K)';
    Xdi=Xd(:,K);
    for pp=1:length(Pd2)
        J=find(~isnan(Xdi(pp,:)));
        if isempty(J), continue, end
        if min(abs(Ld(K(J))-Ld2(nt)))<1.5*corr,
            Edi=exp(- ( (Ldi(pp,:)-Ld2(nt))/corr ).^2 );
            Xd2(pp,nt)=nansum(Xdi(pp,:).*Edi,2)./nansum((Xdi(pp,:)*0+1).*Edi,2);
        end
    end
end
zt2=Xd2;

% lissage fluo
Ld=time;Pd=pres;Xd=zsi;
Ld2=(1:corr:floor(max(Ld)))';Pd2=pres;Xd2=NaN*zeros(length(Pd2),length(Ld2));
for nt=1:length(Ld2),
    K=find(abs(Ld-Ld2(nt))<5*corr);
    if isempty(K), continue, end
    Ldi=ones(size(Pd))*Ld(K)';
    Xdi=Xd(:,K);
    for pp=1:length(Pd2)
        J=find(~isnan(Xdi(pp,:)));
        if isempty(J), continue, end
        if min(abs(Ld(K(J))-Ld2(nt)))<1.5*corr,
            Edi=exp(- ( (Ldi(pp,:)-Ld2(nt))/corr ).^2 );
            Xd2(pp,nt)=nansum(Xdi(pp,:).*Edi,2)./nansum((Xdi(pp,:)*0+1).*Edi,2);
        end
    end
end
zs2=Xd2;

[t,p]=meshgrid(Ld2,Pd2);

% trace transects
if exist('hfig','var') & hfig~=0
    figure(hfig); clf
elseif hfig==0
    hfig = figure('visible','off');
else
    hfig = figure;
end
colormap(jet)

H=[];

h2=subplot(2,1,1);hold on,xlabel('days'),ylabel('depth')
h4=subplot(2,1,2);hold on,xlabel('days'),ylabel('depth')
pos=get(gca,'position');
axes('position',[pos(1) .96 pos(3) .02]); axis off;
pcolor([Ccontour;Ccontour]); shading flat
set(gca,'fontsize',8,'ytick',[],'xtick',[]);

set(gcf,'CurrentAxes',h2);cla,hold on,box on
if any(~isnan(zt2(:))) & length(Ld2)>1
    h22=pcolor(t,p,zt2);shading flat;H=[H h22];%colorbar;
    [Cs,h22]=contour(t,p,zt2,Tcontour,'linecolor','k');H=[H h22];
    set(gca,'XLim',[min(Ccontour) max(Ccontour)],'Ylim',[0 pmax],...
        'fontsize',8,'TickDir','out','ydir','reverse');
    caxis([min(Tcontour) max(Tcontour)]);
    %clabel(Cs,h22,'manual', 'HorizontalAlignment','center','BackgroundColor','w');
    %clabel(Cs,h22,'manual');
    %set(h2,'XTickLabel',datestr(get(gca,'XTick')+693962,'dd mmm'));
    clabel(Cs,h22,'fontsize',8,'labelspacing',500);
    H=[H plot(time,Pbot,'.k')];
    xlim=get(gca,'xlim');ylim=get(gca,'ylim');
    text(xlim(1)+.05*(xlim(2)-xlim(1)),ylim(1)+.95*(ylim(2)-ylim(1)),'IN-SITU TEMP');
else
    set(gcf,'CurrentAxes',h2), axis off
end

if any(~isnan(zs2(:))) & length(Ld2)>1
    set(gcf,'CurrentAxes',h4),cla,hold on,box on
    h23=pcolor(t,p,zs2);shading flat;H=[H h23];%colorbar;
    [Cs,h23]=contour(t,p,zs2,Fcontour,'linecolor','k');H=[H h23];
    set(gca,'XLim',[min(Ccontour) max(Ccontour)],'Ylim',[0 pmax],...
        'fontsize',8,'TickDir','out','ydir','reverse');
    caxis([min(Fcontour) max(Fcontour)]);
    %clabel(Cs,h23,'manual');
    clabel(Cs,h23,'fontsize',8,'labelspacing',500);
    %set(h4,'XTickLabel',datestr(get(gca,'XTick')+693962,'dd mmm'));
    H=[H plot(time,Pbot,'.k')];
    xlim=get(gca,'xlim');ylim=get(gca,'ylim');
    text(xlim(1)+.05*(xlim(2)-xlim(1)),ylim(1)+.95*(ylim(2)-ylim(1)),'CHLA');
else
    set(gcf,'CurrentAxes',h4), axis off
end
    
% print figure
if exist('printname','var') & length(printname)>4 & strcmp(printname(end-3),'.')
    format_figure_centred([28 21]);
    suffix=printname(end-2:end);
    if strcmp(suffix,'png')
        feval(@print,'-dpng','-r300',printname);
    elseif strcmp(suffix,'eps')
        feval(@print,'-depsc2','-tiff',printname);
    else
        feval(@print,['-d' suffix],'-r300',printname);
    end
elseif exist('printname','var')
    error('can''t read printing format');
end

close(hfig);
