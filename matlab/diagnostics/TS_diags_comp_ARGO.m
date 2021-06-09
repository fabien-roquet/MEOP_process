function H=TSdiags_comp_ARGO(argo_qc,argo_qc_raw,argo_clim,conf)
% function TSdiags_ARGO
%   plot TS plot diagnostics for a set of profiles
%
% INPUT:
% argo_qc: struct variable obtained using the function ARGO_load_qc
% conf: structure with configuration parameters
%   Z: vector of size (np,1) used to color lines.
%   edges: edges used in function histc to determine colors. (if NaN, unique color)
%   str: descriptive string used in the caption
%   lim=[south, north, west, east] used to set limits of spatial map.
%   Tlim=[minT maxT] used to set limits of potential temperature maps.
%   Slim=[minS maxS] used to set limits of salinity maps.
%   hfig: figure handle number. if 0, a new figure is opened.
%   nomfig: print name. no print if empty string. suffix indicates the print
%       format.
%
% OUTPUT
% H: vector list of graphical handles

if argo_qc.np==0
    return
end
warning off

H=[];

if ~isempty(argo_clim),
    P_clim=argo_clim.PRES_ADJUSTED;
    T_clim=argo_clim.TEMP_ADJUSTED;
    S_clim=argo_clim.PSAL_ADJUSTED;
    Np_clim=argo_clim.nr+1;

end

P_raw=argo_qc_raw.PRES;
T_raw=argo_qc_raw.TEMP;
S_raw=argo_qc_raw.PSAL;

P=argo_qc.PRES;
T=argo_qc.TEMP;
S=argo_qc.PSAL;
argo_qc.Tmask=double(argo_qc.TEMP_QC<2); argo_qc.TEMP(argo_qc.TEMP_QC>2)=NaN;
argo_qc.Smask=double(argo_qc.PSAL_QC<2); argo_qc.PSAL(argo_qc.PSAL_QC>2)=NaN;

N=argo_qc.np;
Np=argo_qc.nr+1;

if ~isfield(conf,'Z')
    conf.Z=argo_qc.index_tag;
    conf.edges=1:max(conf.Z);
end

if ~isfield(conf,'str'),
    conf.str ={};
end

Ncol=length(conf.edges);

if Ncol<2, % no coloring case
    n=length(N);
    bin = ones(N,1);
    col=[0 0 1];
else % general case
    [n,bin] = histc(conf.Z,conf.edges);
    col=jet(Ncol);
end

if ~isfield(conf,'Tlim') | length(conf.Tlim)~=2
    conf.Tlim=[-2 6];
end
if ~isfield(conf,'Slim') | length(conf.Slim)~=2
    conf.Slim=[33.5 35];
end

if isfield(conf,'hfig') & conf.hfig~=0
    figure(conf.hfig); clf
else
    figure;conf.hfig=gcf;
end

if ~isfield(conf,'list_tag'),
    conf.list_tag = unique(bin);
end
Iqc=find(ismember(argo_qc.index_tag,conf.list_tag));

if ~isfield(conf,'lim')
    
    lon=argo_qc.LONGITUDE(Iqc); is_lon_center_180=0;
    if any(lon<45&lon>-45),
        mlon=floor(min(lon)/10)*10; Mlon=ceil(max(lon)/10)*10;
    else
        lon(lon<0)=lon(lon<0)+360;
        mlon=floor(min(lon)/5)*5-5; Mlon=ceil(max(lon)/5)*5+5;
        is_lon_center_180=1;
    end
    conf.lim=[floor(min(argo_qc.LATITUDE(Iqc))/2)*2-2 ceil(max(argo_qc.LATITUDE(Iqc))/2)*2+2 mlon Mlon];
end

%% initialisation

ax1=axes('position',[0 .55 .25 .4]);axis off

ax2=axes('position',[.25 .55 .3 .4]);
set(gca,'fontsize',7);
xlabel('THETA-S RAW PROFILES','fontsize',8);
hold on
set_tsdiag(conf.Slim,conf.Tlim,0,27:.2:28);
if ~isempty(argo_clim),
    % plot clim TS data
    J=find( ...
        argo_clim.LATITUDE>=conf.lim(1) & argo_clim.LATITUDE<=conf.lim(2) & ...
        argo_clim.LONGITUDE>=conf.lim(3) & argo_clim.LONGITUDE<=conf.lim(4) ...
        );
    positions = [argo_clim.LONGITUDE(J),argo_clim.LATITUDE(J)] ;
    [positions,ia,ic] = unique(positions,'rows'); J=J(ia);
    DT = delaunayTriangulation(positions);
    vi = nearestNeighbor(DT,[argo_qc.LONGITUDE(Iqc),argo_qc.LATITUDE(Iqc)]);
    I=find(ismember(DT.Points(vi,:),positions,'rows'));
    J=J(unique(vi));
    nn=length(J);
    if nn
        Pts=[P_clim(:,J);zeros(1,nn)*NaN]; Pts=Pts(:);
        Tts=[T_clim(:,J);zeros(1,nn)*NaN];  Tts=Tts(:);
        Sts=[S_clim(:,J);zeros(1,nn)*NaN];  Sts=Sts(:);
        It=find(~isnan(Tts));It=union(It,Np_clim:Np_clim:length(Tts));
        Is=find(~isnan(Sts));Is=union(Is,Np_clim:Np_clim:length(Sts));
        Pt=Pts(It);Tt=Tts(It);Ps=Pts(Is);Ts=Tts(Is);Ss=Sts(Is);
        Ts_pot=sw_ptmp(Ss,Ts,Ps,0);
        plot(Ss,Ts_pot,'-','color','k');
    end
end
for ss=[34.6:.05:34.8]
    plot([ss ss],[conf.Tlim(1) conf.Tlim(2)],'k')
end

ax3=axes('position',[.65 .55 .3 .4]);
set(gca,'fontsize',7);
xlabel('THETA-S COR PROFILES','fontsize',8);
hold on
set_tsdiag(conf.Slim,conf.Tlim,0,27:.2:28);
if ~isempty(argo_clim) & nn,
    plot(Ss,Ts_pot,'-','color','k');
end
for ss=[34.6:.05:34.8]
    plot([ss ss],[conf.Tlim(1) conf.Tlim(2)],'k')
end

ax4=axes('position',[.1 .05 .5 .42]);
axis(conf.lim([3:4 1:2]));
hold on
if ~isempty(argo_clim) & nn,
    plot(argo_clim.LONGITUDE(J),argo_clim.LATITUDE(J),'.',...
        'color','k','markersize',10)
end

descr_str = conf.str;
base_descr_str = descr_str;
H2=[];

for jj=conf.list_tag,
    
    J=find(vector(bin)==jj & ...
        argo_qc.LATITUDE>=conf.lim(1) & argo_qc.LATITUDE<=conf.lim(2) & ...
        argo_qc.LONGITUDE>=conf.lim(3) & argo_qc.LONGITUDE<=conf.lim(4) ...
        );
    nn=length(J);
    
    if nn,

        % plot raw TS data
        Pts=[P_raw(:,J);zeros(1,nn)*NaN]; Pts=Pts(:);
        Tts=[T_raw(:,J);zeros(1,nn)*NaN];  Tts=Tts(:);
        Sts=[S_raw(:,J);zeros(1,nn)*NaN];  Sts=Sts(:);
        It=find(~isnan(Tts));It=union(It,Np:Np:length(Tts));
        Is=find(~isnan(Sts));Is=union(Is,Np:Np:length(Sts));
        Pt=Pts(It);Tt=Tts(It);Ps=Pts(Is);Ts=Tts(Is);Ss=Sts(Is);
        Ts_pot=sw_ptmp(Ss,Ts,Ps,0);
        
        axes(ax2);
        H=[H plot(Ss,Ts_pot,'-','color',col(jj,:))];

        % plot cor TS data
        Pts=[P(:,J);zeros(1,nn)*NaN]; Pts=Pts(:);
        Tts=[T(:,J);zeros(1,nn)*NaN];  Tts=Tts(:);
        Sts=[S(:,J);zeros(1,nn)*NaN];  Sts=Sts(:);
        It=find(~isnan(Tts));It=union(It,Np:Np:length(Tts));
        Is=find(~isnan(Sts));Is=union(Is,Np:Np:length(Sts));
        Pt=Pts(It);Tt=Tts(It);Ps=Pts(Is);Ts=Tts(Is);Ss=Sts(Is);
        Ts_pot=sw_ptmp(Ss,Ts,Ps,0);
        
        axes(ax3);
        H=[H plot(Ss,Ts_pot,'-','color',col(jj,:))];
        
        axes(ax4);
        H=[H plot(argo_qc.LONGITUDE(J),argo_qc.LATITUDE(J),'.',...
            'color',col(jj,:),'markersize',10)];
        
        axes(ax1);
        descr_str={descr_str{:},' ', ...
            sprintf('%s: %d profiles',argo_qc.platform_number{J(1)},nn)};
        
        if isfield(conf,'pause') & conf.pause
            pause
        end
        
    end
    
end

H2=text(0.5,0.5,descr_str, ...
    'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);


% print figure
if isfield(conf,'nomfig') & length(conf.nomfig)>4 & strcmp(conf.nomfig(end-3),'.')
    format_figure_centred([28 21]);
    suffix=conf.nomfig(end-2:end);
    if strcmp(suffix,'png')
        feval(@print,'-dpng','-r300',conf.nomfig);
    elseif strcmp(suffix,'eps')
        feval(@print,'-depsc2','-tiff',conf.nomfig);
    else
        feval(@print,['-d' suffix],'-r300',conf.nomfig);
    end
    H=[];
end



