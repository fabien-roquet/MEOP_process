function[h,hc]=topoplot(varargin)
%TOPOPLOT  Plot regional or global topography at one-sixth degree resolution.
%   __________________________________________________________________
%
%   *|* topoplot.png --- Filled contour plot of the continents and shelves.  
%   Type 'jhelp topoplot' to view this image. *|*
%   __________________________________________________________________
%
%   TOPOPLOT with no input arguments makes a patch plot of the earth's
%   surface topography within the axis of the current plot.  Continents are
%   black and the continental shelf to a depth of 500~m is gray.
%
%   TOPOPLOT CONTINENTS plots the continents only as a patch plots with
%   the continents in black.  
%
%   If the current plot is empty, TOPOPLOT uses the domain -180 to 180 
%   degrees in longitude, and -80.666 to 80.666 degrees in latitude.
%
%   These two options use PATCHCONTOURF, which does not use the current 
%   colormap and so allow one to allow plot a colored field in addition to
%   the topography.  The remaining options use CONTOURF.
%   __________________________________________________________________
% 
%   TOPOPLOT(REGION,V) makes a filled contour plot of the earth's surface
%   topography within a region specified by REGION, at contour intervals V.   
%
%   REGION is an array with the format [WEST EAST SOUTH NORTH]. Longitudes 
%   may either be specified on the interval [-180, 180] or on [0, 360].
%
%   V is an array of contour levels, in kilometers, with positive values 
%   for heights above sea level and negative values for below sea level.
%
%   By default, TOPOPLOT uses a flipped grayscale colormap, with high
%   elevations colored dark gray, and low elevations colored light gray.
%   TOPOPLOT uses CONTOURF with the contours themselves not shown. 
%  
%   TOPOPLOT uses the one-sixth degree resolution topography file JTOPO.MAT
%   included with JLAB, based on the Smith and Sandwell dataset.  The 
%   latitude range +/- 80.666 is covered. Type 'help jtopo' for details. 
%
%   TOPOPLOT(REGION) with V unspecified uses V=[0 -1/2 -1], showing the
%   continents in black and the contintental shelf at 500 m depth in gray.
%
%   To plot the continents in black, use TOPOPLOT(REGION,[-1/2 0]).
%
%   TOPOPLOT([],V) sets REGION to the axis limits of the current plot.
%   __________________________________________________________________
%
%   Additional options
%
%   TOPOPLOT(REGION,V,VC) additionally draws the individual contour levels
%   in VC, with the color shading still based on the array V.
%
%   TOPOPLOT(REGION,[],VC) draws only the contours in VC, with no shading.
%   By default, these will be drawn as heavy white lines.
%  
%   TOPOPLOT(REGION,V,VC,STY) uses linestyle STY for the individual contour
%   levels.  STY is in LINESTYLE form, e.g. STY='2g--'.
%
%   [H,HC]=TOPOPLOT(REGION,V,VC) returns the handles to the contours
%   corresponding to levels V and VC, respectively.
%   
%   TOPOPLOT(...,'M_MAP') will work with Rich Pawlowicz's M_MAP package by
%   calling M_CONTOUR and M_CONTOURF.
%   __________________________________________________________________
%
%   See also JTOPO, REGIONPLOT, READTOPO.
% 
%   'topoplot --f' generates the above sample figures and some others.
%
%   Usage: topoplot(region,v);
%          topoplot(region,v,'m_map');
%          [h,hc]=topoplot(region,v,vc);
%          [h,hc]=topoplot(region,v,vc,'2m--');
%          [h,hc]=topoplot(region,v,vc,'2m--','m_map');
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2015 J.M. Lilly --- type 'help jlab_license' for details
 

%
%   TOPOPLOT with no input arguments sets REGION to the axis limits of the
%   current plot using the default value of V=[0 -1/2 -1]. 
%
%   With no input arguments, if the current axes is empty, the axis is 
%   first set to [-180 180 -80.666 80.666], topography is plotted with the
%   default choice of V, and HOLD ON is applied.  Try 'FIGURE, TOPOPLOT'. 

if nargin>0
    if strcmpi(varargin{1}, '--f')
       makefigs_topoplot,return
    end
end


sty='2w';
str='matlab';
optionstr='   ';

for i=1:3
    if length(varargin)>=1
        if ischar(varargin{end})
            if strfind(lower(varargin{end}),'m_m')
                str=lower(varargin{end});
            elseif ~isempty(strfind(lower(varargin{end}),'con'))||~isempty(strfind(lower(varargin{end}),'she'))
                optionstr=varargin{end};
            else
                sty=varargin{end};
            end
            varargin=varargin(1:end-1);
        end
    end
end
ax=[];
region=[];

if length(varargin)>0
    region=varargin{1};
end

if isempty(get(gcf,'children'))
   axis([-180 180 -80.666 80.666])
   hold on
   bnewplot=true;
else
   bnewplot=false;
end

if isempty(region)
    region=axis;
end
ax=region;

if ~bnewplot&&((ax(2)-ax(1))<360-1/3)
    region(1)=region(1)-1/6;
    region(2)=region(2)+1/6;
    region(3)=region(3)-1/6;
    region(4)=region(4)+1/6;
end

depths=[];
if length(varargin)>1
    depths=varargin{2};
end

depthsc=[];
if length(varargin)>2
    if ~ischar(varargin{3})
        depthsc=varargin{3};
        varargin=varargin(4:end);
    else
        varargin=varargin(3:end);
    end
end

if length(depths)==1
    depths=[depths depths];
end
if length(depthsc)==1
    depthsc=[depthsc depthsc];
end

load jtopo
use jtopo

west=region(1);
east=region(2);
south=region(3);
north=region(4);

a=find(lat>south,1,'first');
b=find(lat<north,1,'last');
if a>1,a=a-1;end  %Expand the box by one point N/S, for plotting purposes.
if b<length(lat),b=b+1;end 

lat=lat(a:b);
topo=topo(a:b,:);

[lon,topo]=lonshift(west,lon,topo);

%Repeat last column
lon=[lon lon(end)+lon(end)-lon(end-1)];
topo=[topo topo(:,1)];

a=find(lon>west,1,'first');
b=find(lon>=east,1,'first');

if a>1,a=a-1;end %Expand the box by one point E/W, for plotting purposes.
if b<length(lon),b=b+1;end
        
lon=lon(a:b);
topo=topo(:,a:b);

hold on
h=[];
hc=[];

if strcmpi(optionstr(1:3),'she')
    patchcontourf(lon,lat,topo,-2,[1 1 1]*2/3,str);
    hold on
    patchcontourf(lon,lat,topo,0,'k',str);
elseif strcmpi(optionstr(1:3),'con')
    patchcontourf(lon,lat,topo,0,'k',str);
else
    if strcmpi(str(1:3),'m_m')
        if ~isempty(depths)
            for kdepth=1:length(depths)
                patchcontourf(lon,lat,topo,depths(kdepth),'k',str);
                hold on
            end
        end
        if ~isempty(depthsc)
            if length(depthsc)==1
                depthsc=[depthsc depthsc];
            end
            [c,hc]=m_contour(lon,lat,topo,-depthsc,'linecolor','k');
            hh=get(hc,'children');
            linestyle(hh,sty);
        end
        if isempty(depths) & isempty(depthsc)
            patchcontourf(lon,lat,topo,-2,[1 1 1]*2/3,str);
            hold on
            patchcontourf(lon,lat,topo,0,'k',str);
        end
    else
        if ~isempty(depths)
            for kdepth=1:length(depths)
                patchcontourf(lon,lat,topo,depths(kdepth),'k',str);
                hold on
            end
        end
        if ~isempty(depthsc)
            if length(depthsc)==1
                depthsc=[depthsc depthsc];
            end
            [c,hc]=contour(lon,lat,topo,-depthsc,'linecolor','k');
            hh=get(hc,'children');
            linestyle(hh,sty);
        end
        if isempty(depths) & isempty(depthsc)
            patchcontourf(lon,lat,topo,-2,[1 1 1]*2/3,str);
            hold on
            patchcontourf(lon,lat,topo,0,'k',str);
        end
        if ~isempty(ax)
            axis(ax);
        end
    end
    if ~isempty(depths), caxis([min(depths) max(depths)]); end
end


%If I have used contourf, drop the z-levels
if ~isempty(h)
    hh = get(h,'Children');   
    for i=1:length(hh)
        zdata = ones(size( get(hh(i),'XData') ));
        set(hh(i), 'ZData',-10*zdata)
    end
end

boxon

if isempty(get(get(gca,'xlabel'),'string'))
    xlabel('Longitude')
end
if isempty(get(get(gca,'ylabel'),'string'))
    ylabel('Latitude')
end
%figure,pcolor(lon,lat,topo),shading interp

if nargout == 0
    clear h hc
end

