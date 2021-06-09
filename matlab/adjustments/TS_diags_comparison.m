function conf = TSdiags_comparison(conf)
% function TSdiags_ARGO
%   plot TS plot diagnostics for a set of profiles
%
% INPUT:
%
% conf: structure with configuration parameters
%   argo_qc: struct variable obtained using the function ARGO_load_qc
%   argo_wod : struct variable (WOD data) (optional)
%   argo_cora : struct variable (CORA data) (optional)
%   argo_meop: struct variable (additional qc-ed meop data) (optional)
%   argo_qc2: struct variable (other tags from same deployment) (optional)
%   Z: vector of size (np,1) used to color lines.
%   edges: edges used in function histc to determine colors. (if NaN, unique color)
%   str: descriptive string used in the caption
%   lim=[south, north, west, east] used to set limits of spatial map.
%   Tlim=[minT maxT] used to set limits of potential temperature maps.
%   Slim=[minS maxS] used to set limits of salinity maps.
%   hfig: figure handle number. if 0, the figure is created but not visible.
%   nomfig: print name. no print if empty string. suffix indicates the print
%       format.
%
% OUTPUT
% H: vector list of graphical handles

if isfield(conf,'hfig') & conf.hfig~=0
    figure(conf.hfig); clf
elseif conf.hfig==0
    conf.hfig = figure('visible','off');
else
    conf.hfig = figure;
end


if conf.argo_qc.np==0
    close(conf.hfig)
    return
end
warning off

if ~isfield(conf,'list_tag'),
    conf.list_tag = conf.argo_qc.list_descr;
end

Iqc=find(ismember(conf.argo_qc.platform_number,conf.list_tag) & sum(conf.argo_qc.PSAL_QC==1)'>1);
if isempty(Iqc),
    %disp(['no valid profile: ' conf.nomfig]);
    close(conf.hfig)
    return
end

if ~isfield(conf,'lim')
    lon=conf.argo_qc.LONGITUDE(Iqc);
    lon(lon>180)=lon(lon>180)-360;
    zone1=any(lon>-180&lon<-90);
    zone2=any(lon> -90&lon<  0);
    zone3=any(lon>   0&lon< 90);
    zone4=any(lon>  90&lon<180);
    if ~ (zone1 & zone4) | (zone2 & zone3),
        mlon=floor(min(lon)/10)*10; Mlon=ceil(max(lon)/10)*10; is_lon_center_180=0;
    else
        lon(lon<0)=lon(lon<0)+360;
        mlon=floor(min(lon)/10)*10; Mlon=ceil(max(lon)/10)*10; is_lon_center_180=1;
    end
    conf.lim=[floor(min(conf.argo_qc.LATITUDE(Iqc))/2)*2-2 ceil(max(conf.argo_qc.LATITUDE(Iqc))/2)*2+2 mlon Mlon];
else
    conf.lim(conf.lim>180)=conf.lim(conf.lim>180)-360;
end



H=[];

argo_wod =0;
data_wod.strdata='wod';
if isfield(conf,'argo_wod') & ~isempty(conf.argo_wod)
    argo_wod=1;
    data_wod.LATITUDE=conf.argo_wod.LATITUDE;
    data_wod.LONGITUDE=conf.argo_wod.LONGITUDE;
    data_wod.P=conf.argo_wod.PRES;
    data_wod.T=conf.argo_wod.TEMP;
    data_wod.S=conf.argo_wod.PSAL;
    data_wod.Np=conf.argo_wod.nr+1;
    
    data_wod.col='k';
end

argo_cora =0;
data_cora.strdata='cora';
if isfield(conf,'argo_cora') & ~isempty(conf.argo_cora)
    argo_cora=1;
    data_cora.LATITUDE=conf.argo_cora.LATITUDE;
    data_cora.LONGITUDE=conf.argo_cora.LONGITUDE;
    data_cora.P=conf.argo_cora.PRES;
    data_cora.T=conf.argo_cora.TEMP;
    data_cora.S=conf.argo_cora.PSAL;
    data_cora.Np=conf.argo_cora.nr+1;
    
    data_cora.col='k';
end

argo_meop=0;
data_meop.strdata='meop';
if isfield(conf,'argo_meop') & ~isempty(conf.argo_meop)
    argo_meop=1;
    data_meop.LATITUDE=conf.argo_meop.LATITUDE;
    data_meop.LONGITUDE=conf.argo_meop.LONGITUDE;
    data_meop.P=conf.argo_meop.PRES;
    data_meop.T=conf.argo_meop.TEMP;
    data_meop.S=conf.argo_meop.PSAL;
    data_meop.Np=conf.argo_meop.nr+1;
    data_meop.col='r';
    data_meop.PLATFORM_NUMBER=conf.argo_meop.PLATFORM_NUMBER;
end

argo_qc2=0;
data_qc2.strdata='tags';
if isfield(conf,'argo_qc2') & ~isempty(conf.argo_qc2)
    argo_qc2=1;
    data_qc2.LATITUDE=conf.argo_qc2.LATITUDE;
    data_qc2.LONGITUDE=conf.argo_qc2.LONGITUDE;
    data_qc2.P=conf.argo_qc2.PRES;
    data_qc2.T=conf.argo_qc2.TEMP;
    data_qc2.S=conf.argo_qc2.PSAL;
    data_qc2.Np=conf.argo_qc2.nr+1;
    data_qc2.col='c';
    data_qc2.PLATFORM_NUMBER=conf.argo_qc2.PLATFORM_NUMBER;
    if ~isfield(conf,'nomfig2'),
        argo_qc2=0;
        disp('name file conf_adjustment.nomfig2 missing. No figure created');
    end
end

if isfield(conf,'noclim') & conf.noclim
    argo_wod =0;
    argo_cora =0;
    argo_meop=0;
end


N=conf.argo_qc.np;
Ntag=conf.argo_qc.ntag;
Np=conf.argo_qc.nr+1;

P = conf.argo_qc.PRES;
T = conf.argo_qc.TEMP;
S = conf.argo_qc.PSAL;
T(conf.argo_qc.TEMP_QC>2)=NaN;
S(conf.argo_qc.PSAL_QC>2)=NaN;

conf.isoffset = 'No variable offset';
P_raw = P;
T_raw = T*NaN;
S_raw = S*NaN;
if Ntag == 1
    T_raw = T + conf.T1 * 1e-3 * P_raw + conf.T2;
    S_raw = S + conf.S1 * 1e-3 * P_raw + conf.S2;
    if conf.offset
        S_raw = S_raw - repmat(conf.offset,size(S,1),1);
        conf.isoffset = sprintf('Variable offset %g +/- %g',mean(conf.offset),std(conf.offset));
    end
end

    
if ~isfield(conf,'Z')
    conf.Z=conf.argo_qc.index_tag;
    conf.edges=1:max(conf.Z);
end

if ~isfield(conf,'str'),
    conf.str ={};
end

Ncol=length(conf.edges);

if Ncol<2, % no coloring case
    n=length(Ntag);
    bin = ones(Ntag,1);
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


%% initialisation

ax1=axes('position',[0 .55 .25 .4]);axis off

ax2=axes('position',[.05 .55 .4 .4]);
set(gca,'fontsize',7);
xlabel('THETA-S RAW PROFILES','fontsize',8);
hold on
set_tsdiag(conf.Slim,conf.Tlim,0,25:.2:28);
for ss=[round(conf.Slim(1),1):.1:round(conf.Slim(2),1)]
    plot([ss ss],[conf.Tlim(1) conf.Tlim(2)],'g','linewidth',.3)
end
for ss=[34.6:.05:34.8]
    plot([ss ss],[conf.Tlim(1) conf.Tlim(2)],'m','linewidth',.3)
end
for tt=[round(conf.Tlim(1)):round(conf.Tlim(2))]
    plot([conf.Slim(1) conf.Slim(2)],[tt tt],'g','linewidth',.3)
end

if argo_wod,
    data=data_wod;
    % plot clim TS data
    if conf.lim(3)<conf.lim(4)
        J=find( ...
            data.LATITUDE>=conf.lim(1) & data.LATITUDE<=conf.lim(2) & ...
            data.LONGITUDE>=conf.lim(3) & data.LONGITUDE<=conf.lim(4) ...
            );
    else
        aux = data.LONGITUDE;
        aux(aux<0)=aux(aux<0)+360;
        J=find( ...
            data.LATITUDE>=conf.lim(1) & data.LATITUDE<=conf.lim(2) & ...
            aux>=conf.lim(4) & aux<=conf.lim(3) ...
            );
    end
    if length(J)==0
        argo_wod=0;
    else
        positions = [data.LONGITUDE(J),data.LATITUDE(J)] ;
        [positions,ia,ic] = unique(positions,'rows'); J=J(ia);
        if exist('delaunayTriangulation')
            DT = delaunayTriangulation(positions);
        else
            DT = DelaunayTri(positions);
        end
        if ~isempty(DT.ConnectivityList)
            vi = nearestNeighbor(DT,[conf.argo_qc.LONGITUDE(Iqc),conf.argo_qc.LATITUDE(Iqc)]);
            data.J=J(unique(vi));
            data.nn=length(data.J);
        else
            data.nn=0;
            data.J=[];
        end
        if data.nn
            nn=data.nn; J=data.J;
            Pts=[data.P(:,J);zeros(1,nn)*NaN];  Pts=Pts(:);
            Tts=[data.T(:,J);zeros(1,nn)*NaN];  Tts=Tts(:);
            Sts=[data.S(:,J);zeros(1,nn)*NaN];  Sts=Sts(:);
            It=find(~isnan(Tts));It=union(It,data.Np:data.Np:length(Tts));
            Is=find(~isnan(Sts));Is=union(Is,data.Np:data.Np:length(Sts));
            data.Pt=Pts(It);
            data.Tt=Tts(It);
            data.Ps=Pts(Is);
            data.Ts=Tts(Is);
            data.Ss=Sts(Is);
            data.Ts_pot=sw_ptmp(data.Ss,data.Ts,data.Ps,0);
            plot(data.Ss,data.Ts_pot,'-','color',data.col,'linewidth',.3);
        else
            argo_wod=0;
        end
        data_wod=data;
    end
end

if argo_cora,
    data=data_cora;
    % plot clim TS data
    if conf.lim(3)<conf.lim(4)
        J=find( ...
            data.LATITUDE>=conf.lim(1) & data.LATITUDE<=conf.lim(2) & ...
            data.LONGITUDE>=conf.lim(3) & data.LONGITUDE<=conf.lim(4) ...
            );
    else
        aux = data.LONGITUDE;
        aux(aux<0)=aux(aux<0)+360;
        J=find( ...
            data.LATITUDE>=conf.lim(1) & data.LATITUDE<=conf.lim(2) & ...
            aux>=conf.lim(4) & aux<=conf.lim(3) ...
            );
    end
    if length(J)==0
        argo_cora=0;
    else
        positions = [data.LONGITUDE(J),data.LATITUDE(J)] ;
        [positions,ia,ic] = unique(positions,'rows'); J=J(ia);
        if exist('delaunayTriangulation')
            DT = delaunayTriangulation(positions);
        else
            DT = DelaunayTri(positions);
        end
        if ~isempty(DT.ConnectivityList)
            vi = nearestNeighbor(DT,[conf.argo_qc.LONGITUDE(Iqc),conf.argo_qc.LATITUDE(Iqc)]);
            data.J=J(unique(vi));
            data.nn=length(data.J);
        else
            data.nn=0;
            data.J=[];
        end
        if data.nn
            nn=data.nn; J=data.J;
            Pts=[data.P(:,J);zeros(1,nn)*NaN];  Pts=Pts(:);
            Tts=[data.T(:,J);zeros(1,nn)*NaN];  Tts=Tts(:);
            Sts=[data.S(:,J);zeros(1,nn)*NaN];  Sts=Sts(:);
            It=find(~isnan(Tts));It=union(It,data.Np:data.Np:length(Tts));
            Is=find(~isnan(Sts));Is=union(Is,data.Np:data.Np:length(Sts));
            data.Pt=Pts(It);
            data.Tt=Tts(It);
            data.Ps=Pts(Is);
            data.Ts=Tts(Is);
            data.Ss=Sts(Is);
            data.Ts_pot=sw_ptmp(data.Ss,data.Ts,data.Ps,0);
            plot(data.Ss,data.Ts_pot,'-','color',data.col,'linewidth',.3);
        else
            argo_cora=0;
        end
        data_cora=data;
    end
end

if argo_meop,
    data=data_meop;
    % plot clim TS data
    if conf.lim(3)<conf.lim(4)
        J=find( ...
            data.LATITUDE>=conf.lim(1) & data.LATITUDE<=conf.lim(2) & ...
            data.LONGITUDE>=conf.lim(3) & data.LONGITUDE<=conf.lim(4) & ...
            ~ismember(cellstr(data.PLATFORM_NUMBER'),...
            unique(cellstr(conf.argo_qc.PLATFORM_NUMBER(:,Iqc)'))) ...
            );
    else
        aux = data.LONGITUDE;
        aux(aux<0)=aux(aux<0)+360;
        J=find( ...
            data.LATITUDE>=conf.lim(1) & data.LATITUDE<=conf.lim(2) & ...
            aux>=conf.lim(4) & aux<=conf.lim(3) & ...
            ~ismember(cellstr(data.PLATFORM_NUMBER'),...
            unique(cellstr(conf.argo_qc.PLATFORM_NUMBER(:,Iqc)'))) ...
            );
    end
    if length(J)==0
        argo_meop=0;
    else
        positions = [data.LONGITUDE(J),data.LATITUDE(J)] ;
        [positions,ia,ic] = unique(positions,'rows'); J=J(ia);
        if exist('delaunayTriangulation')
            DT = delaunayTriangulation(positions);
        else
            DT = DelaunayTri(positions);
        end
        if ~isempty(DT.ConnectivityList)
            vi = nearestNeighbor(DT,[conf.argo_qc.LONGITUDE(Iqc),conf.argo_qc.LATITUDE(Iqc)]);
            data.J=J(unique(vi));
            data.nn=length(data.J);
        else
            data.nn=0;
            data.J=[];
        end
        if data.nn
            nn=data.nn; J=data.J;
            Pts=[data.P(:,J);zeros(1,nn)*NaN];  Pts=Pts(:);
            Tts=[data.T(:,J);zeros(1,nn)*NaN];  Tts=Tts(:);
            Sts=[data.S(:,J);zeros(1,nn)*NaN];  Sts=Sts(:);
            It=find(~isnan(Tts));It=union(It,data.Np:data.Np:length(Tts));
            Is=find(~isnan(Sts));Is=union(Is,data.Np:data.Np:length(Sts));
            data.Pt=Pts(It);
            data.Tt=Tts(It);
            data.Ps=Pts(Is);
            data.Ts=Tts(Is);
            data.Ss=Sts(Is);
            data.Ts_pot=sw_ptmp(data.Ss,data.Ts,data.Ps,0);
            plot(data.Ss,data.Ts_pot,'-','color',data.col,'linewidth',.3);
        else
            argo_meop=0;
        end
        data_meop=data;
    end
end

if argo_wod==0 & argo_cora==0 & argo_meop==0
    %disp(['no reference profiles: ' conf.nomfig]);
    close(conf.hfig)
    return
end

ax3=axes('position',[.55 .55 .4 .4]);
set(gca,'fontsize',7);
xlabel('THETA-S COR PROFILES','fontsize',8);
hold on
set_tsdiag(conf.Slim,conf.Tlim,0,25:.2:28);
for ss=[round(conf.Slim(1),1):.1:round(conf.Slim(2),1)]
    plot([ss ss],[conf.Tlim(1) conf.Tlim(2)],'g','linewidth',.3)
end
for ss=[34.6:.05:34.8]
    plot([ss ss],[conf.Tlim(1) conf.Tlim(2)],'m','linewidth',.3)
end
for tt=[round(conf.Tlim(1)):round(conf.Tlim(2))]
    plot([conf.Slim(1) conf.Slim(2)],[tt tt],'g','linewidth',.3)
end
if argo_wod & data_wod.nn,
    data=data_wod;
    plot(data.Ss,data.Ts_pot,'-','color',data.col,'linewidth',.3);
end
if argo_cora & data_cora.nn,
    data=data_cora;
    plot(data.Ss,data.Ts_pot,'-','color',data.col,'linewidth',.3);
end
if argo_meop & data_meop.nn,
    data=data_meop;
    plot(data.Ss,data.Ts_pot,'-','color',data.col,'linewidth',.3);
end


ax4=axes('position',[.05 .05 .4 .42]);
axis(conf.lim([3:4 1:2]));
hold on
if argo_wod & data_wod.nn,
    data=data_wod;
    plot(data.LONGITUDE(data.J),data.LATITUDE(data.J),'.',...
        'color',data.col,'markersize',10)
end
if argo_cora & data_cora.nn,
    data=data_cora;
    plot(data.LONGITUDE(data.J),data.LATITUDE(data.J),'.',...
        'color',data.col,'markersize',10)
end
if argo_meop & data_meop.nn,
    data=data_meop;
    plot(data.LONGITUDE(data.J),data.LATITUDE(data.J),'.',...
        'color',data.col,'markersize',10)
end

descr_str = conf.str;
base_descr_str = descr_str;
H2=[];

if isfield(conf,'pause') & conf.pause
    pause
end

for jj=1:length(conf.list_tag),
    
    if conf.lim(3)<conf.lim(4)
        J=find(bin==jj & sum(conf.argo_qc.PSAL_QC==1)'>1 & ...
            conf.argo_qc.LATITUDE>=conf.lim(1) & conf.argo_qc.LATITUDE<=conf.lim(2) & ...
            conf.argo_qc.LONGITUDE>=conf.lim(3) & conf.argo_qc.LONGITUDE<=conf.lim(4) ...
            );
    else
        aux = conf.argo_qc.LONGITUDE;
        aux(aux<0)=aux(aux<0)+360;
        J=find(bin==jj & sum(conf.argo_qc.PSAL_QC==1)'>1 & ...
            conf.argo_qc.LATITUDE>=conf.lim(1) & conf.argo_qc.LATITUDE<=conf.lim(2) & ...
            aux>=conf.lim(4) & aux<=conf.lim(3) ...
            );
    end
    nn=length(J);
    
    if nn,
        
        % plot raw TS data
        set(gcf,'CurrentAxes',ax2);
        for ii=1:length(J),
            try
                I1 = find(~isnan(P_raw(:,J(ii))));
                I2 = find(~isnan(P_raw(:,J(ii)).*T_raw(:,J(ii))));
                T_raw(I1,J(ii)) = interp1(P_raw(I2,J(ii)),T_raw(I2,J(ii)),P_raw(I1,J(ii)));
                I3 = find(~isnan(P_raw(:,J(ii)).*S_raw(:,J(ii))));
                S_raw(I1,J(ii)) = interp1(P_raw(I3,J(ii)),S_raw(I3,J(ii)),P_raw(I1,J(ii)));
            end
        end
        Pts=[P_raw(:,J);zeros(1,nn)*NaN]; Pts=Pts(:);
        Tts=[T_raw(:,J);zeros(1,nn)*NaN];  Tts=Tts(:);
        Sts=[S_raw(:,J);zeros(1,nn)*NaN];  Sts=Sts(:);
        Its=find(~isnan(Tts)|~isnan(Sts));Its=union(Its,Np:Np:length(Pts));
        It=find(~isnan(Tts));It=union(It,Np:Np:length(Tts));
        Is=find(~isnan(Sts));Is=union(Is,Np:Np:length(Sts));
        Pt=Pts(It);Tt=Tts(It);Ps=Pts(Is);Ts=Tts(Is);Ss=Sts(Is);
        Ts_pot=sw_ptmp(Ss,Ts,Ps,0);
        H=[H plot(Ss,Ts_pot,'-','color',col(jj,:),'linewidth',.3)];
        
        % plot cor TS data
        set(gcf,'CurrentAxes',ax3);
        for ii=1:length(J),
            try
                I1 = find(~isnan(P(:,J(ii))));
                I2 = find(~isnan(P(:,J(ii)).*T(:,J(ii))));
                T(I1,J(ii)) = interp1(P(I2,J(ii)),T(I2,J(ii)),P(I1,J(ii)));
                I3 = find(~isnan(P(:,J(ii)).*S(:,J(ii))));
                S(I1,J(ii)) = interp1(P(I3,J(ii)),S(I3,J(ii)),P(I1,J(ii)));
            end
        end
        Pts=[P(:,J);zeros(1,nn)*NaN]; Pts=Pts(:);
        Tts=[T(:,J);zeros(1,nn)*NaN];  Tts=Tts(:);
        Sts=[S(:,J);zeros(1,nn)*NaN];  Sts=Sts(:);
        It=find(~isnan(Tts));It=union(It,Np:Np:length(Tts));
        Is=find(~isnan(Sts));Is=union(Is,Np:Np:length(Sts));
        Pt=Pts(It);Tt=Tts(It);Ps=Pts(Is);Ts=Tts(Is);Ss=Sts(Is);
        Ts_pot=sw_ptmp(Ss,Ts,Ps,0);
        H=[H plot(Ss,Ts_pot,'-','color',col(jj,:),'linewidth',.3)];
        
        set(gcf,'CurrentAxes',ax4);
        H=[H plot(conf.argo_qc.LONGITUDE(J),conf.argo_qc.LATITUDE(J),'.',...
            'color',col(jj,:),'markersize',10)];
        
        ax6=axes('position',[.55 .05 .4 .42]);axis off
        set(gcf,'CurrentAxes',ax6);
        if isfield(conf.argo_qc,'list_smru_platform')
            descr_str={descr_str{:},' ', ...
                sprintf('%s :\n %d profiles \n T1 = %g ; T2 =%g \n S1 = %g ; S2=%g \n %s \n time : %s \n %s',...
                conf.argo_qc.list_smru_platform{jj},nn,conf.T1(jj),conf.T2(jj),conf.S1(jj),conf.S2(jj)),conf.isoffset,...
                datestr(conf.argo_qc.JULD(1)),datestr(conf.argo_qc.JULD(end))};
        else
            descr_str={descr_str{:},' ', ...
                sprintf('%s: %d profiles',conf.argo_qc.platform_number{J(1)},nn)};
        end
        
        if isfield(conf,'pause') & conf.pause
            pause
        end
        
    end
    
end

H2=text(0.5,0.5,descr_str, ...
    'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);

% print figure
if isfield(conf,'nomfig') & length(conf.nomfig)>4 & strcmp(conf.nomfig(end-3),'.')
    format_figure_centred([20 15]);
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


%% comparison figure with tags from same deployment
if ~argo_qc2,
    close(conf.hfig)
    return;
end

clf

ax1=axes('position',[0 .55 .25 .4]);axis off

ax2=axes('position',[.05 .55 .4 .4]);
set(gca,'fontsize',7);
xlabel('THETA-S RAW PROFILES','fontsize',8);
hold on
set_tsdiag(conf.Slim,conf.Tlim,0,27:.2:28);
for ss=[round(conf.Slim(1),1):.1:round(conf.Slim(2),1)]
    plot([ss ss],[conf.Tlim(1) conf.Tlim(2)],'g','linewidth',.3)
end
for ss=[34.6:.05:34.8]
    plot([ss ss],[conf.Tlim(1) conf.Tlim(2)],'m','linewidth',.3)
end
for tt=[round(conf.Tlim(1)):round(conf.Tlim(2))]
    plot([conf.Slim(1) conf.Slim(2)],[tt tt],'g','linewidth',.3)
end

data=data_qc2;
% plot clim TS data
if conf.lim(3)<conf.lim(4)
    J=find( ...
        data.LATITUDE>=conf.lim(1) & data.LATITUDE<=conf.lim(2) & ...
        data.LONGITUDE>=conf.lim(3) & data.LONGITUDE<=conf.lim(4) & ...
        ~ismember(cellstr(data.PLATFORM_NUMBER'),...
        unique(cellstr(conf.argo_qc.PLATFORM_NUMBER(:,Iqc)'))) ...
        );
else
    aux = data.LONGITUDE;
    aux(aux<0)=aux(aux<0)+360;
    J=find( ...
        data.LATITUDE>=conf.lim(1) & data.LATITUDE<=conf.lim(2) & ...
        aux>=conf.lim(4) & aux<=conf.lim(3) & ...
        ~ismember(cellstr(data.PLATFORM_NUMBER'),...
        unique(cellstr(conf.argo_qc.PLATFORM_NUMBER(:,Iqc)'))) ...
        );
end
if length(J)==0
    argo_qc2=0;
else
    positions = [data.LONGITUDE(J),data.LATITUDE(J)] ;
    [positions,ia,ic] = unique(positions,'rows'); J=J(ia);
    if exist('delaunayTriangulation')
        DT = delaunayTriangulation(positions);
    else
        DT = DelaunayTri(positions);
    end
    if ~isempty(DT.ConnectivityList)
        vi = nearestNeighbor(DT,[conf.argo_qc.LONGITUDE(Iqc),conf.argo_qc.LATITUDE(Iqc)]);
        data.J=J(unique(vi));
        data.nn=length(data.J);
    else
        data.nn=0;
        data.J=[];
    end
    if data.nn
        nn=data.nn; J=data.J;
        for ii=1:length(J),
            try
                I1 = find(~isnan(data.P(:,J(ii))));
                I2 = find(~isnan(data.P(:,J(ii)).*data.T(:,J(ii))));
                data.T(I1,J(ii)) = interp1(data.P(I2,J(ii)),data.T(I2,J(ii)),data.P(I1,J(ii)));
                I3 = find(~isnan(data.P(:,J(ii)).*data.S(:,J(ii))));
                data.S(I1,J(ii)) = interp1(data.P(I3,J(ii)),data.S(I3,J(ii)),data.P(I1,J(ii)));
            end
        end
        Pts=[data.P(:,J);zeros(1,nn)*NaN];  Pts=Pts(:);
        Tts=[data.T(:,J);zeros(1,nn)*NaN];  Tts=Tts(:);
        Sts=[data.S(:,J);zeros(1,nn)*NaN];  Sts=Sts(:);
        It=find(~isnan(Tts));It=union(It,data.Np:data.Np:length(Tts));
        Is=find(~isnan(Sts));Is=union(Is,data.Np:data.Np:length(Sts));
        data.Pt=Pts(It);
        data.Tt=Tts(It);
        data.Ps=Pts(Is);
        data.Ts=Tts(Is);
        data.Ss=Sts(Is);
        data.Ts_pot=sw_ptmp(data.Ss,data.Ts,data.Ps,0);
        plot(data.Ss,data.Ts_pot,'-','color',data.col,'linewidth',.3);
    else
        argo_qc2=0;
    end
    data_qc2=data;
end

if ~argo_qc2,
    close(conf.hfig)
    return;
end

ax3=axes('position',[.55 .55 .4 .4]);
set(gca,'fontsize',7);
xlabel('THETA-S COR PROFILES','fontsize',8);
hold on
set_tsdiag(conf.Slim,conf.Tlim,0,27:.2:28);
for ss=[round(conf.Slim(1),1):.1:round(conf.Slim(2),1)]
    plot([ss ss],[conf.Tlim(1) conf.Tlim(2)],'g','linewidth',.3)
end
for ss=[34.6:.05:34.8]
    plot([ss ss],[conf.Tlim(1) conf.Tlim(2)],'m','linewidth',.3)
end
for tt=[round(conf.Tlim(1)):round(conf.Tlim(2))]
    plot([conf.Slim(1) conf.Slim(2)],[tt tt],'g','linewidth',.3)
end
data=data_qc2;
plot(data.Ss,data.Ts_pot,'-','color',data.col,'linewidth',.3);


ax4=axes('position',[.05 .05 .4 .42]);
axis(conf.lim([3:4 1:2]));
hold on
data=data_qc2;
plot(data.LONGITUDE(data.J),data.LATITUDE(data.J),'.',...
    'color',data.col,'markersize',10)

descr_str = conf.str;
base_descr_str = descr_str;
H2=[];

if isfield(conf,'pause') & conf.pause
    pause
end

for jj=1:length(conf.list_tag),
    
    if conf.lim(3)<conf.lim(4)
        J=find(bin==jj & sum(conf.argo_qc.PSAL_QC==1)'>1 & ...
            conf.argo_qc.LATITUDE>=conf.lim(1) & conf.argo_qc.LATITUDE<=conf.lim(2) & ...
            conf.argo_qc.LONGITUDE>=conf.lim(3) & conf.argo_qc.LONGITUDE<=conf.lim(4) ...
            );
    else
        aux = conf.argo_qc.LONGITUDE;
        aux(aux<0)=aux(aux<0)+360;
        J=find(bin==jj & sum(conf.argo_qc.PSAL_QC==1)'>1 & ...
            conf.argo_qc.LATITUDE>=conf.lim(1) & conf.argo_qc.LATITUDE<=conf.lim(2) & ...
            aux>=conf.lim(4) & aux<=conf.lim(3) ...
            );
    end
    nn=length(J);
    
    if nn,
        
        % plot raw TS data
        set(gcf,'CurrentAxes',ax2);
        for ii=1:length(J),
            try
                I1 = find(~isnan(P_raw(:,J(ii))));
                I2 = find(~isnan(P_raw(:,J(ii)).*T_raw(:,J(ii))));
                T_raw(I1,J(ii)) = interp1(P_raw(I2,J(ii)),T_raw(I2,J(ii)),P_raw(I1,J(ii)));
                I3 = find(~isnan(P_raw(:,J(ii)).*S_raw(:,J(ii))));
                S_raw(I1,J(ii)) = interp1(P_raw(I3,J(ii)),S_raw(I3,J(ii)),P_raw(I1,J(ii)));
            end
        end
        Pts=[P_raw(:,J);zeros(1,nn)*NaN]; Pts=Pts(:);
        Tts=[T_raw(:,J);zeros(1,nn)*NaN];  Tts=Tts(:);
        Sts=[S_raw(:,J);zeros(1,nn)*NaN];  Sts=Sts(:);
        Its=find(~isnan(Tts)|~isnan(Sts));Its=union(Its,Np:Np:length(Pts));
        It=find(~isnan(Tts));It=union(It,Np:Np:length(Tts));
        Is=find(~isnan(Sts));Is=union(Is,Np:Np:length(Sts));
        Pt=Pts(It);Tt=Tts(It);Ps=Pts(Is);Ts=Tts(Is);Ss=Sts(Is);
        Ts_pot=sw_ptmp(Ss,Ts,Ps,0);
        H=[H plot(Ss,Ts_pot,'-','color',col(jj,:),'linewidth',.3)];
        
        % plot cor TS data
        set(gcf,'CurrentAxes',ax3);
        for ii=1:length(J),
            try
                I1 = find(~isnan(P(:,J(ii))));
                I2 = find(~isnan(P(:,J(ii)).*T(:,J(ii))));
                T(I1,J(ii)) = interp1(P(I2,J(ii)),T(I2,J(ii)),P(I1,J(ii)));
                I3 = find(~isnan(P(:,J(ii)).*S(:,J(ii))));
                S(I1,J(ii)) = interp1(P(I3,J(ii)),S(I3,J(ii)),P(I1,J(ii)));
            end
        end
        Pts=[P(:,J);zeros(1,nn)*NaN]; Pts=Pts(:);
        Tts=[T(:,J);zeros(1,nn)*NaN];  Tts=Tts(:);
        Sts=[S(:,J);zeros(1,nn)*NaN];  Sts=Sts(:);
        It=find(~isnan(Tts));It=union(It,Np:Np:length(Tts));
        Is=find(~isnan(Sts));Is=union(Is,Np:Np:length(Sts));
        Pt=Pts(It);Tt=Tts(It);Ps=Pts(Is);Ts=Tts(Is);Ss=Sts(Is);
        Ts_pot=sw_ptmp(Ss,Ts,Ps,0);
        H=[H plot(Ss,Ts_pot,'-','color',col(jj,:),'linewidth',.3)];
        
        set(gcf,'CurrentAxes',ax4);
        H=[H plot(conf.argo_qc.LONGITUDE(J),conf.argo_qc.LATITUDE(J),'.',...
            'color',col(jj,:),'markersize',10)];
        
        ax6=axes('position',[.55 .05 .4 .42]);axis off
        if isfield(conf.argo_qc,'list_smru_platform')
            descr_str={descr_str{:},' ', ...
                sprintf('%s :\n %d profiles \n T1 = %g ; T2 =%g \n S1 = %g ; S2=%g \n %s \n time : %s \n %s',...
                conf.argo_qc.list_smru_platform{jj},nn,conf.T1(jj),conf.T2(jj),conf.S1(jj),conf.S2(jj)),conf.isoffset,...
                datestr(conf.argo_qc.JULD(1)),datestr(conf.argo_qc.JULD(end))};
        else
            descr_str={descr_str{:},' ', ...
                sprintf('%s: %d profiles',conf.argo_qc.platform_number{J(1)},nn)};
        end
        
        if isfield(conf,'pause') & conf.pause
            pause
        end
        
    end
    
end

H2=text(0.5,0.5,descr_str, ...
    'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',12);

% print figure
if isfield(conf,'nomfig2') & length(conf.nomfig2)>4 & strcmp(conf.nomfig2(end-3),'.')
    format_figure_centred([20 15]);
    suffix=conf.nomfig2(end-2:end);
    if strcmp(suffix,'png')
        feval(@print,'-dpng','-r300',conf.nomfig2);
    elseif strcmp(suffix,'eps')
        feval(@print,'-depsc2','-tiff',conf.nomfig2);
    else
        feval(@print,['-d' suffix],'-r300',conf.nomfig2);
    end
    H=[];
end

close(conf.hfig);


