function histogram_ARGO(argo_qc,hfig,printname)
% function histogram_ARGO
%   plot histograms diagnosing for a set of profiles
%
% INPUT:
% argo_qc: argo_qc in a struct variable, obtained through the function ARGO_load_qc
% hfig: figure handle number. if 0, the figure is not visible.
% printname: print name. no print if empty string. suffix indicates the print
%       format.


if argo_qc.np==0
    return
end
warning off

P=argo_qc.PRES;
T=argo_qc.TEMP;
S=argo_qc.PSAL;

D=sw_dens(S,T,0)-1000;
N=argo_qc.np;

if exist('hfig','var') & hfig~=0
    figure(hfig);clf
elseif hfig == 0
    hfig = figure('visible','off');
else
    hfig=figure;
end

Pbot=argo_qc.JULD_LOCATION*NaN;
Pmask=double(~isnan(P));
Tmask=double(~isnan(T));
Smask=double(~isnan(S));
Ibot=sum(Pmask,1);
for kk=1:N, 
    if Ibot(kk), Pbot(kk)=P(Ibot(kk),kk); end, 
end
year=argo_qc.JULD_LOCATION(:)/365.25;

% initialisation
ax1=subplot(2,3,1);
[n,xout]=hist(Pbot(sum(Smask.*Tmask,1)~=0),0:200:1000);
[n2,xout2]=hist(Pbot(sum(Smask+Tmask,1)~=0),0:200:1000);
bar(xout',[n2;n-n2]','stacked');set(gca,'xlim',[0 1000]);
title('BOTTOM PRES.');

ax2=subplot(2,3,2);
[n,xout]=hist(year(sum(Smask.*Tmask,1)~=0),2003.5:1:2022.5);
[n2,xout2]=hist(year(sum(Smask+Tmask,1)~=0),2003.5:1:2022.5);
bar(xout',[n2;n-n2]','stacked');
set(gca,'xlim',[2004 2022],'xtick',2004.5:1:2022.5, ...
    'xticklabel',{'','05','','','','','10','','','','','15','','','','','20','',''})
title('DATE')

ax3=subplot(2,3,3);
month=(year-floor(year))*12;
[n,xout]=hist(month(sum(Smask.*Tmask,1)~=0),.5:12.5);
[n2,xout2]=hist(month(sum(Smask+Tmask,1)~=0),.5:12.5);
bar(xout',[n2;n-n2]','stacked');
set(gca,'xlim',[0 12],'xtick',.5:11.5, ...
    'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'})
title('MONTH')

ax4=subplot(2,3,4);
hist(P(~isnan(T)),0:200:1000);set(gca,'xlim',[0 1000])
title('PRESSURE')

ax5=subplot(2,3,5);
hist(T(~isnan(T)),-2:1:30);set(gca,'xlim',[-3 30])
title('IN-SITU TEMP')

ax6=subplot(2,3,6);
hist(S(~isnan(S)),0:1:40);set(gca,'xlim',[0 40])
title('SALI')

suptitle(sprintf('%d T+S profiles, %d T and/or S profiles',...
    length(find(sum(Smask.*Tmask,1)~=0)), length(find(sum(Smask+Tmask,1)~=0))));

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

close(hfig)

