function[h]=regionplot(varargin)
%REGIONPLOT  Plots a simple box indicating a latitude / longitude region.
%
%   REGIONPLOT(REGION) draws a latitude / longitue box indicated by REGION
%   on the current axes. 
% 
%   REGION is an array with the format [WEST EAST SOUTH NORTH]. Longitudes 
%   may either be specified on the interval [-180, 180] or on [0, 360].
%
%   By default, longitudes will be converted to [-180,180] unless the
%   region overlaps the dateline, in which case [0,360] will be used.
%
%   The current axes will be held as the box is plotted, then returned to 
%   its previous HOLD state. 
%
%   REGION can also be a cell array of N arrays, each with the format 
%   REGION{N}=[WEST EAST SOUTH NORTH].  All N regions will then be plotted.
%
%   REGIONPLOT is compatible with the M_MAP toolbox, as described below.
%
%   H=REGIONPLOT outputs the line handle H to the plotted box.
%   __________________________________________________________________
%
%   Additional options
%
%   H=REGIONPLOT returns the graphics handle to the plotted box. 
%
%   REGIONPLOT(REGION,STY) uses style string STY, in the LINESTYLE format,
%   with a default value of STY='k'.
%
%   REGIONPLOT(REGION,'M_MAP') will work with Rich Pawlowicz's M_MAP 
%   package by calling M_PLOT.  
% 
%   The additional input parameters can be input together in either order,
%   i.e. REGIONPLOT(REGION,STY,'M_MAP') or REGIONPLOT(REGION,'M_MAP',STY).
%   __________________________________________________________________
%
%   See also INREGION, FLOATREGION, TOPOPLOT.
%
%   'regionplot --f' generates some sample figures.
%
%   Usage: regionplot(region);
%          h=regionplot(region,'m_map'); 
%          h=regionplot(region,'m_map','2k'); 
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2014 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--f')
    regionplot_figure,return
end
region=varargin{1};
varargin=varargin(2:end);
sty='k';
str='matlab';
for i=1:2
    if length(varargin)>=1
        if ischar(varargin{1})
            if strfind(lower(varargin{1}),'m_m')
                str=lower(varargin{1});
            else
                sty=varargin{1};
            end
            varargin=varargin(2:end);
        end
    end
end

if ~iscell(region)
    h=regionplot_one(region,str,sty);
    set(h,'visible','on')
    %set(h(~isnan(h)),'visible','on')
else
    bhold=ishold;
    set(gca,'visible','off')
    hold on
    storestate=get(gcf,'BackingStore');
    set(gcf,'BackingStore','off')
    for i=1:length(region)
         h(i)=regionplot_one(region{i},str,sty);
    end
    linestyle(h,'default');
    set(h,'visible','on')
    %set(h(~isnan(h)),'visible','on')
    
    set(gcf,'BackingStore',storestate)
    if ~bhold
        hold off
    end
end

if nargout==0
    clear h
end


function[h]=regionplot_one(region,str,sty)
            
west=region(1);
east=region(2);
south=region(3);
north=region(4);

ax=axis;
if ax(1)>180 || deg180(west)>deg180(east)
    west=deg360(west);
    east=deg360(east);
else 
    west=deg180(west);
    east=deg180(east);
end

bhold=ishold;
hold on
if strcmpi(str(1:3),'m_m')
     h=m_plot([west east east west west],[south south north north south],'visible','off');
else
     h=plot([west east east west west],[south south north north south],'visible','off');
end
linestyle(h,sty);
if ~bhold,
    hold off
end


function[]=regionplot_figure
 
load ebasnfloats
use ebasnfloats
region=[-30 -21 24 35];
figure,cellplot(lon,lat),latratio,regionplot(region),axis tight
title('Floats from Eastern Basin experiment')
xlabel('Longitude'),ylabel('Latitude'),boxon
ax=axis;
fontsize 14 12 12 12

if exist('m_map')==7
    figure,
    m_proj('albers equal-area conic','lon',ax(1:2),'lat',ax(3:4));
    cellplot(lon,lat,'m_map')
    m_grid('linestyle','none')
    regionplot(region,'m_map')
    title('Floats from Eastern Basin experiment using M-MAP')
    fontsize 14 12 12 12
else
    disp('REGIONPLOT skipping M_MAP example figure since M_MAP not detected.')
end

%use ebasnfloats
%%Region shifted to be greater than 180
%region=[-30+250 -21+250 24 35];
%%Region shifted to overlap dateline
%region=[-30+205 -21+205 24 35];
%figure,cellplot(celladd(205,lon),lat),latratio,regionplot(region),axis tight
%title('Floats from Eastern Basin experiment')
%xlabel('Longitude'),ylabel('Latitude'),boxon
%ax=axis;


%reporttest('REGIONPLOT',aresame())
