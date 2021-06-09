function H=TSdiags_ARGO(argo_qc,Z,edges,descr_str,maplim,Tlim,Slim,Dlim,hfig,printname)
% function TSdiags_ARGO
%   plot TS plot diagnostics for a set of profiles
%
% INPUT:
% argo_qc: struct variable obtained using the function ARGO_load_qc
% Z: vector of size (np,1) used to color lines.
% edges: edges used in function histc to determine colors. (if NaN, unique color)
% descr_str: descriptive string used in the caption
% maplim=[south, north, west, east] used to set limits of spatial map.
% Tlim=[minT maxT] used to set limits of potential temperature maps.
% Slim=[minS maxS] used to set limits of salinity maps.
% Dlim=[minD maxD] used to set limits of sigma0 maps.
% hfig: figure handle number. if 0, a new figure is opened.
% printname: print name. no print if empty string. suffix indicates the print
%       format.
%
% OUTPUT
% H: vector list of graphical handles


if argo_qc.np==0
    return
end
warning off

H=[];

P=argo_qc.PRES;
T=argo_qc.TEMP;
isSal=0;if isfield(argo_qc,'PSAL'), isSal=1; S=argo_qc.PSAL; end
isfluo=0; if isfield(argo_qc,'CHLA'), isfluo=1; F=argo_qc.CHLA; end
    
D=sw_dens(S,T,0)-1000;
N=argo_qc.np;
Np=argo_qc.nr+1;
Ncol=length(edges);

[n,bin] = histc(Z,edges);
col=jet(Ncol);

if ~exist('Tlim','var') | length(Tlim)~=2
    Tlim=[min(-2,floor(min(T(:)))) max(10,ceil(max(T(:))))];
end
if ~exist('Slim','var') | length(Slim)~=2
    Slim=[min(32.5,floor(min(S(:)))) max(35,ceil(max(S(:))))];
end
if ~exist('Dlim','var') | length(Dlim)~=2
    Dlim=[min(26,floor(min(D(:)))) max(28.5,ceil(max(D(:))))];
end
if ~exist('Flim','var') | length(Dlim)~=2
    Flim=[0 2];
end

if exist('hfig','var') & hfig~=0
    figure(hfig);clf
elseif hfig == 0
    hfig = figure('visible','off');
else
    hfig=figure;
end

%% initialisation
if ~strcmp(get(hfig,'Tag'),'fiche_init')
    clf,set(hfig,'Tag','fiche_init')
    if maplim(3)<maplim(4), mloc=maplim(3):20:maplim(4);
    else mloc=[fliplr(-(-180:20:-maplim(3))) -180:20:maplim(4)];
    end
    
    ax1=axes('position',[0 .55 .25 .4]);axis off
    ax2=axes('position',[.25 .55 .2 .4]);set(gca,'YDir','reverse','fontsize',7);xlabel('IN-SITU TEMP. PROFILES','fontsize',8);hold on;grid on
    ax3=axes('position',[.5 .55 .2 .4]);set(gca,'YDir','reverse','fontsize',7);xlabel('SALI. PROFILES','fontsize',8);hold on;grid on
    ax4=axes('position',[.75 .55 .2 .4]);set(gca,'YDir','reverse','fontsize',7);xlabel('SIGMA0 PROFILES','fontsize',8);hold on;grid on
    ax5=axes('position',[.1 .05 .5 .42]);axis off
    m_proj('Equidistant Cylindrical',...
        'longitude',maplim(3:4),'latitude',maplim(1:2));
    m_grid('box','on','xtick',mloc,'xticklabels',mloc,...
        'ytick',maplim(1):10:maplim(2),'yticklabels',maplim(1):10:maplim(2),...
        'tickdir','out','linestyle','none',...
        'XaxisLocation','bottom','YaxisLocation','left');
    m_gshhs_l('patch',[.7 .7 .7],'linestyle','none');
    hold on
    ax6=axes('position',[.65 .05 .3 .4]);set(gca,'fontsize',7);xlabel('THETA-S PROFILES','fontsize',8);hold on;grid on;
    
else
    ax1=findobj('position',[0 .55 .25 .4]);
    ax2=findobj('position',[.25 .55 .2 .4]);
    ax3=findobj('position',[.5 .55 .2 .4]);
    ax4=findobj('position',[.75 .55 .2 .4]);
    ax5=findobj('position',[.1 .05 .5 .42]);
    ax6=findobj('position',[.65 .05 .3 .4]);
end

set(gcf,'CurrentAxes',ax1);
H=[H text(0.5,0.5,descr_str, ...
    'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12)];


for jj=1:length(n), % bug if only one color !
    J=find(bin==jj);
    if n(jj),
        
        for ii=1:length(J),
            try
            I1 = find(~isnan(P(:,J(ii))));
            I2 = find(~isnan(P(:,J(ii)).*T(:,J(ii))));
            T(I1,J(ii)) = interp1(P(I2,J(ii)),T(I2,J(ii)),P(I1,J(ii)));
            I3 = find(~isnan(P(:,J(ii)).*S(:,J(ii))));
            S(I1,J(ii)) = interp1(P(I3,J(ii)),S(I3,J(ii)),P(I1,J(ii)));
            end
        end
        Pts=[P(:,J);zeros(1,n(jj))*NaN];  Pts=Pts(:);
        Tts=[T(:,J);zeros(1,n(jj))*NaN];  Tts=Tts(:);
        Sts=[S(:,J);zeros(1,n(jj))*NaN];  Sts=Sts(:);
        It=find(~isnan(Tts));It=union(It,Np:Np:length(Tts));
        Is=find(~isnan(Sts));Is=union(Is,Np:Np:length(Sts));
        Pt=Pts(It);Tt=Tts(It);Ps=Pts(Is);Ts=Tts(Is);Ss=Sts(Is);
        Ts_pot=sw_ptmp(Ss,Ts,Ps,0);
        Ds=sw_dens(Ss,Ts_pot,0)-1000;
        
        set(gcf,'CurrentAxes',ax2);H=[H plot(Tt,Pt,'-','color',col(jj,:))];set(gca,'xlim',Tlim)
        set(gcf,'CurrentAxes',ax3);H=[H plot(Ss,Ps,'-','color',col(jj,:))];set(gca,'xlim',Slim)
        set(gcf,'CurrentAxes',ax4);H=[H plot(Ds,Ps,'-','color',col(jj,:))];set(gca,'xlim',Dlim)
        set(gcf,'CurrentAxes',ax5);
        H=[H m_plot(argo_qc.LONGITUDE(J),argo_qc.LATITUDE(J),'.',...
                    'color',col(jj,:),'markersize',5)];
        
        set(gcf,'CurrentAxes',ax6);
        H=[H plot(Ss,Ts_pot,'-','color',col(jj,:))];
                
    end
end

set(gcf,'CurrentAxes',ax6);H=[H set_tsdiag(Slim,Tlim,0)];


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
    delete(H);H=[];
elseif exist('printname','var')
    error('can''t read printing format');
end



