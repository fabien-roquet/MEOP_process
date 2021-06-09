function H=tsdiagram(S,T,P,sigma);
% TSDIAGRM  This function plots a temp vs. salinity diagram,
%            with selected density contours.
%
%            TSDIAGRM(S,T,P,SIGMA) draws contours lines of density anomaly
%            SIGMA (kg/m^3) at pressure P (dbars), given a range of
%            salinity (ppt) and temperature (deg C) in the 2-element vectors
%            S,T.  The freezing point (if visible) will be indicated.
%
%            TSDIAGRM(S,T,P) draws several randomly chosen contours.


%Adapted from RP (WHOI) 9/Dec/91
%                 7/Nov/92 Changed for Matlab 4.0
%                 14/Mar/94 Made P optional.

if nargin==1,
    error('tsdiagram: Not enough calling parameters');
elseif nargin==2,
    P=0;
elseif nargin==0,
    S=get(gca,'XLim');
    T=get(gca,'YLim');
    P=0;
end;

Sg=S(1)+[0:50]/50*(S(2)-S(1));
Tg=T(1)+[0:50]'/50*(T(2)-T(1));
SG=sw_dens(ones(size(Tg))*Sg,Tg*ones(size(Sg)),P(1))-1000;
if nargin~=4
    Z=SG; mZ=min(SG(:)); MZ=max(SG(:)); sigma=list_contours(mZ,MZ,5);
end
axis([S(1) S(2) T(1) T(2)]);
[CS,h]=contour(Sg,Tg,SG,sigma,'Linecolor','k','Linewidth',.5);

%plot freezing temp.
freezeT=sw_fp(S,P(1));
h2=plot(S,freezeT,'k--','Linewidth',.5);

%% Label with pressure, then return to other axes
% set(gca,'FontSize',9);
% xlabel('Salinity','FontSize',9);
% ylabel('Potential Temperature ^oC','FontSize',9);
% text(S(1),T(2), ...
% [' Pressure = ' int2str(P(1)) ' dbars'],'horiz','left','Vert','top');
%clabel(CS,h,'manual','fontsize',7);
clabel(CS,h,'fontsize',7);


H=[h h2];

