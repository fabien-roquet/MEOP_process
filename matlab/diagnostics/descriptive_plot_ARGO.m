function descriptive_plot_ARGO(argo_data,hfig,str_title,printname,conf)
% function descriptive_plot_ARGO
%   plot profile positions + histograms diagnosing for a set of profiles
%
% INPUT:
% argo_data: argo_data in a struct variable, obtained through the function
%   ARGO_load
% hfig: figure handle number. if 0, a new figure is opened.
% str_title: descriptive string used in the title
% printname: print name. no print if empty string. suffix indicates the print
%       format.


if argo_data.np==0
    return
end
warning off



% data
lat=argo_data.LATITUDE;
lon=argo_data.LONGITUDE;
P=argo_data.PRES;
T=argo_data.TEMP;
S=argo_data.PSAL;
D=sw_dens(S,T,0)-1000;
N=argo_data.np;
tag=argo_data.platform_number;
Nprof=length(lon);
Ntag=length(argo_data.list_descr);
GROUP=argo_data.GROUP;
list_group=unique(GROUP);

Pbot=argo_data.JULD*NaN;
Ibot=sum(argo_data.Pmask,1);
for kk=1:N, 
    if Ibot(kk), Pbot(kk)=P(Ibot(kk),kk); end, 
end
year=argo_data.JULD(:)/365.25;
%madt_rio=interpolate_rio09(lat,lon);

% plot
if exist('hfig','var') & hfig~=0
    figure(hfig);
else
    figure;hfig=gcf;
end

ax1=findobj('position',[.13 .11 .475 .75]);
if length(ax1)~=1,
    delete(ax1);
    ax1=axes('position',[.13 .11 .475 .75]);
    if isnumeric(ax1),
        set(gcf,'Tag',num2str(ax1));
    else
        set(gcf,'Tag',num2str(ax1.Position(1)));
    end
    m_proj('gnomonic','lat',-90,'long',0,'radius',50);
    %m_elev('contour',[-2000 -2000],'edgecolor','k');
    m_grid('box','on','xtick',12,...
        'tickdir','out','ytick',[-80:10:-50],...
        'XaxisLocation','top','linest','-');
    m_coast('patch',[.7 .7 .7]);
    hold on
end

%
axes(ax1);
H=[]; Hlegend=[];
col=hsv(length(list_group));
for ii=1:length(list_group),
    Igroup=find(strcmp(GROUP,list_group(ii)));
    list_tag=unique(tag(Igroup));
  %  col=list_color{2}(find(strcmp(conf.list_color(1),list_group(ii))),:);
    for kk=1:length(list_tag),
        I=find(strcmp(list_tag{kk},tag));
        if length(I)>0
            h=m_plot(lon(I),lat(I));
            set(h,'markersize',5,'linewidth',.7,'color',col(ii,:));
            H=[H h];
        end
    end
    Hlegend=[Hlegend h];
end
H2=legend(Hlegend,list_group,'fontsize',10,'location','SW');
set(findobj(get(H2,'Children'),'linewidth',.7),'linewidth',1.5)

ax2=findobj('position',[0 0 .7 1]);
delete(ax2);
ax2=axes('position',[0 0 .7 1],'visible','off');
text(.5,.93,sprintf('(a) %s: %d profiles, %d tags',str_title,Nprof,Ntag),...
    'HorizontalAlignment','center','fontsize',12);
%legend(H,{list_ficseals{Icol,2}},'fontsize',10);

% histo
ax3=axes('position',[.7 .77 .2 .16]);cla,hold on
nlist=[];
for ii=1:length(list_group),
    Igroup=find(strcmp(GROUP,list_group(ii)));
    [n,xout]=hist(Pbot(Igroup),0:100:1500);
    nlist=[nlist;n];
end
if length(list_group)>1, nlist=cumsum(nlist); end
for ii=length(list_group):-1:1,
    b=bar(xout',nlist(ii,:)','stacked');
    set(b,'facecolor',col(ii,:));
end
set(gca,'xlim',[0 1500]);
title('(b) DIVE DEPTH');

ax4=axes('position',[.7 .53 .2 .16]);cla,hold on
nlist=[];
for ii=1:length(list_group),
    Igroup=find(strcmp(GROUP,list_group(ii)));
    [n,xout]=hist(year(Igroup),2003.5:1:2013.5);
    nlist=[nlist;n];
end
if length(list_group)>1, nlist=cumsum(nlist); end
for ii=length(list_group):-1:1,
    b=bar(xout',nlist(ii,:)','stacked');
    set(b,'facecolor',col(ii,:));
end
set(gca,'xlim',[2004 2014],'xtick',2004.5:1:2013.5, ...
    'xticklabel',{'04','05','06','07','08','09','10','11','12','13'})
title('(c) YEAR')

ax5=axes('position',[.7 .29 .2 .16]);cla,hold on
month=(year-floor(year))*12;
nlist=[];
for ii=1:length(list_group),
    Igroup=find(strcmp(GROUP,list_group(ii)));
    [n,xout]=hist(month(Igroup),.5:12.5);
    nlist=[nlist;n];
end
if length(list_group)>1, nlist=cumsum(nlist); end
for ii=length(list_group):-1:1,
    b=bar(xout',nlist(ii,:)','stacked');
    set(b,'facecolor',col(ii,:));
end
set(gca,'xlim',[0 12],'xtick',.5:11.5, ...
    'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'})
title('(d) MONTH')

% data=madt_rio;
% data(madt_rio<-1.1 | isnan(madt_rio))=1;
% data(madt_rio>-1.1 & madt_rio<-.9)=2;
% data(madt_rio>-.9 & madt_rio<-.6)=3;
% data(madt_rio>-.6 & madt_rio<-.1)=4;
% data(madt_rio>-.1)=5;
% ax6=axes('position',[.7 .05 .2 .16]);cla,hold on
% nlist=[];
% for ii=1:length(list_group),
%     Igroup=find(strcmp(GROUP,list_group(ii)));
%     [n,xout]=hist(data(Igroup),1:5);
%     nlist=[nlist;n];
% end
ax6=axes('position',[.7 .04 .2 .16]);cla,hold on
if length(list_group)>1, nlist=cumsum(nlist); end
for ii=length(list_group):-1:1,
    b=bar(xout',nlist(ii,:)','stacked');
    set(b,'facecolor',col(ii,:));
end
set(gca,'xlim',[.5 5.5],'xtick',1:5,'xticklabel',{'SP','SACC','P','SA','ST'});
title('(e) GEO. ZONE')


% ACC front
% axes(ax1);
% ACC_fronts=load_fronts_ACC;
% for kk=2:length(ACC_fronts),
%     ACC_fronts{kk}(ACC_fronts{kk}(:,2)>-40,2)=NaN;
%     h=m_plot(ACC_fronts{kk}(:,1),ACC_fronts{kk}(:,2));
%     set(h,'linewidth',1,'color','k');
%     H=[H h];
% end


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

try, delete(H); end


