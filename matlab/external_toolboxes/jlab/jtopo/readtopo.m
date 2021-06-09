function [mat,lat,lon,offset]=readtopo(varargin)
%READTOPO  Read one-minute topography data from Smith and Sandwell. 
%   _______________________________________________________________________
%   
%   *|* readtopo.png --- Figure of part of the Smith and Sandwell dataset.  
%   Type 'jhelp readtopo' to view this image. *|*
%   _______________________________________________________________________
%
%   [TOPO,LAT,LON]=READTOPO(REGION) extracts regional topography from the
%   the one-minute Smith and Sandwell global database version 18.1.  This 
%   file is included with JDATA, available from at http://www.jmlilly.net.
%
%   REGION is an array of the form REGION=[WEST EAST SOUTH NORTH] which  
%   indicates the *edges* of the desired domain.  Longitudes can be 
%   specified either on the interval [0,360] or the interval [-180,180]. 
%
%   The region may overlap the prime meridian or the dateline.  Region
%   boundaries are interpreted to exclude the poles. 
% 
%   TOPO is a matrix of topography in units of kilometers and is positive 
%   for above sea level, and negative for below sea level.  
%
%   LON is a row vector of longitudes and LAT a column vector of latitudes.
%
%   LAT and LON output by READTOPO are *grid-centered*, that is, they
%   indicate midpoints of the topography cells. LON is uniformly spaced, 
%   but LAT is non-uniformly spaced, with bin sizes decreasing poleward.
%
%   Note that the Smith and Sandwell database is defined only for latitudes 
%   between -80.738 and 80.738.
%   __________________________________________________________________
%
%   Interpolation
%
%   [TOPO,LAT,LON]=READTOPO(REGION,DLAT,DLON) optionally linearly 
%   interpolates the one-minute (1/60 degree) topographic data to a 
%   different resolution, specified by DLAT and DLON, in degrees.   
%
%   READTOPO(REGION,DLAT) with DLON omitted sets DLON=DLAT.
%
%   This is useful for interpolating the high-resolution Smith and Sandwell
%   dataset to a coarser resolution. 
%   __________________________________________________________________
%
%   Alternate path
%
%   By default, READTOPO looks for the topography file topo_18.1.img in 
%   the JLAB directory.  
%
%   [TOPO,LAT,LON]=READTOPO(DIRNAME,REGION...) alternatively specifies  
%   DIRNAME as the path to the parent directory of the file topo_18.1.img. 
%   __________________________________________________________________
%
%   Data and documentation
%
%   The original location for the Smith and Sandwell Global Topography
%    Dataset v. 18.1, file 'topo_18.1.img' is
%
%       http://topex.ucsd.edu/WWW_html/mar_topo.html.
%
%   The reference for the Smith and Sandwell Database is
%
%      Smith, W. H. F., and D. T. Sandwell, Global seafloor topography 
%      from satellite altimetry and ship TOPO soundings, Science, v. 277, 
%      p. 1957-1962, 26 Sept., 1997.
%   __________________________________________________________________
%
%   License and Copyright 
%
%   The data file topo_18.1.img is distributed with JDATA for RESEARCH AND 
%   NON-PROFIT USE ONLY, in accordance with the copyright statement for the
%   Smith and Sandwell dataset.  For details, type 'help topo_copyright'.
%
%   You are free to use and redistribute READTOPO under the terms in the
%   JLAB license, http://www.jmlilly.net/doc/jlab/jlab_license.html.
%
%   READTOPO is self-contained, although tests and sample figure require 
%   JLAB to run.  JLAB is available from http://www.jmlilly.net.
%
%   Send comments, questions, and bug reports to 'eponym@jmlilly.net'.
%   __________________________________________________________________
%    
%   See also JTOPO, REGIONPLOT, and TOPOPLOT in JLAB.
%
%   'readtopo --f' generates the sample figure shown above.
%   'readtopo --t' runs some tests.
%
%   Usage: [topo,lat,lon]=readtopo([west east south north]);
%          [topo,lat,lon]=readtopo(dirname,[west east south north]);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2015 J.M. Lilly --- type 'help jlab_license' for details


%   Note: READTOPO has the same functionality as 'extract_1m' that was 
%   formerly distributed with the Smith and Sandwell databaset, but has 
%   been written from scratch.  Several apparent bugs have been fixed.
%
%   A number of test have been applied to READTOPO.  For example, the calls
%
%        readtopo(dirname,[0 0 80.738 80.738])
%        readtopo(dirname,[360-1/60 360-1/60 -80.738 -80.738]);
%
%   respectively return the first and last points in the file, as expected.

%See http://www.ngdc.noaa.gov/mgg/global/gridregistration.html for a 
%discussion of grid-centered vs. cell-centered.

%From the README at the ftp site... Note Latitude ranges are now +/- 80.738
%and the spacing is now 1/60.  
%
%   The topography is meters above sea level. An even
%   value signifies the cell does not have a ship or coastline
%   measurement while an odd value signifies that it does.
%   The Mercator projected image spans longitudes from 0 E to
%   360 E and latitudes from 72.006 N to -72.006 N.
%   A spherical earth is used for the Mercator projection.
%   The center of the upper left grid cell (i.e. the first
%   integer in the file) is located at 72.0009 N, .01667 E.
%   Longtiudes increase with a 1/30 degree spacing.  The
%   The center of the last integer in the file is located
%   at -72.0009 N, 359.933 E.  Note the latitude spacing is
%   1/30 degree at the equator but decreases as 1/cos(latitude)
%   according to a Mercator projection on a sphere.
%
%Center of first and last latitude bins are at +/- (1/120)/(1/cosd(80.738))
%This matches the tests below. 

if strcmpi(varargin{1}, '--t')
    readtopo_test,return
end
if strcmpi(varargin{1}, '--f')
    readtopo_figure,return
end

if ~ischar(varargin{1})
    dirname=[];
else
    dirname=varargin{1};
    varargin=varargin(2:end);
end

region=varargin{1};

switch length(varargin)
    case 1  
        dlat=0;
        dlon=0;
    case 2 
        dlat=varargin{2};
        dlon=dlat;
    case 3
        dlat=varargin{2};
        dlon=varargin{3};
end
    
minlat=-80.738;
maxlat=80.738;

west=deg180_internal(region(1)-dlon);
east=deg180_internal(region(2)+dlon);

%if east<west,east=east+360;end

deltalon=region(1)-west;

south=max(region(3)-dlat,minlat);
north=min(region(4)+dlat,maxlat);

ddeg=1/60;  
nbytes=21600*2;  % 2-byte integers
skip=17280;

%Some facts: the whole file is 746496000 bytes 
%            there are 17280 different latitude bands so we have
%            746496000/17280=43200 bytes in each stripe of longitude


if north<minlat||south>maxlat
    warning('The Smith and Sandwell database is defined within +/- 80.738 degrees latitude.')
    mat=[];lat=[];lon=[];
else
    if isempty(dirname)
        fid=fopen('topo_18.1.img', 'r','b');  
    else
        fid=fopen([dirname '/topo_18.1.img'], 'r','b');  
    end
    if (fid < 0)
         error(['Could not open database topo_18.1.img.']);
    end

    north=north(:);
    south=south(:);
    east=east(:);
    west=west(:);
    
    if (west<0)&&(east>=0)
        [mat1,lat,lon1,offset1]=readtopo1(fid,south,north,west,360,ddeg,minlat,skip,nbytes);
        [mat2,lat,lon2,offset2]=readtopo1(fid,south,north,0,east,ddeg,minlat,skip,nbytes);
        mat=[mat1 mat2];
        lon=[lon1 lon2];
        offset=[offset1 offset2];
    else
        [mat,lat,lon,offset]=readtopo1(fid,south,north,west,min(east,360-ddeg),ddeg,minlat,skip,nbytes);
    end     
    fclose(fid);

    lon=unwrap((2*pi/360)*lon)*(360/2/pi)+deltalon; %OK.  
    if all(lon(:)>180),
        lon=lon-360;
    end
end
mat=mat/1000;

if nargin>2
    
    %figure,plot(diff(lon)),figure
    %lono=jrad2deg_internal(unwrap(jdeg2rad_internal(lon)));

    [lonog,latog]=meshgrid(lon,lat);
    
    west=deg180_internal(region(1));
    east=deg180_internal(region(2));
    
    if east<west,east=east+360;end
    
    south=region(3);
    north=region(4);

    lat=[south:dlat:north]';
    lon=[west:dlon:east]';
   
    %figure,plot(lon)
    %figure,plot(diff(lono))
    lon1=(2*pi/360)*(unwrap(360/2/pi)*lon);

    [long,latg]=meshgrid(lon,lat);
 
    mat=interp2(lonog,latog,mat,long,latg,'linear');
end


function[mat,lat,lon,offset]=readtopo1(fid,south,north,west,east,ddeg,minlat,skip,nbytes)
[west,east]=deg360_internal(west,east);

arg=log(tand(45+minlat./2)./tand(45+south./2));
%isouth=floor(arg/((2*pi/360)*ddeg)) + skip+1;
%I believe this fixes a problem in EXTRACT_1M
isouth=floor(arg/((2*pi/360)*ddeg)) + skip; 

arg=log(tand(45+minlat./2)./tand(45+north./2));
inorth=floor(arg/((2*pi/360)*ddeg)) + skip+1;
inorth=min(inorth,17280); %Correction for we're right at southern boundary

%In case of only one latitude input
isouth=max(inorth,isouth);
ilat=inorth:isouth;

iwest=ceil(west./ddeg);
iwest=min(iwest,21599);  %Correction for we're right at western boundary

ieast=ceil(east./ddeg)-1;

%In case of only one longitude input
ieast=max(iwest,ieast);
ilon=iwest:ieast;

Nlats=length(ilat);
Nlons=length(ilon);

mat=zeros(Nlats,Nlons);

for i=1:Nlats
        %offset=ilat(i)*nbytes+iwest*2;
        %I believe this fixes another problem in EXTRACT_1M
        offset=(ilat(i)-1)*nbytes+iwest*2;
        fseek(fid, offset, 'bof');
        mat(isouth-ilat(i)+1,:)=fread(fid,[1,Nlons],'int16');
end

lat=zeros(Nlats,1);
argo=tand(45+minlat/2);

arg=exp((2*pi/360).*(ddeg*(skip-(ilat-0.5)))).*argo;
lat(isouth-ilat+1)=2*atand(arg)-90; 
lon=deg180_internal(ddeg*(ilon+0.5)); %To make cell 


function[]=readtopo_figure

if isempty(which('jdata'))
     return
end

% region=[165 -170 54 65];
% [topo1,lat1,lon1]=readtopo(region);
% [topo2,lat2,lon2]=readtopo(region,1/6);
% 
% figure
% subplot(2,2,1)
% pcolor(lon1,lat1,topo1),shading interp, latratio, caxis([-5 2])
% title('Original resolution Smith and Sandwell')
% subplot(2,2,2)
% pcolor(lon2,lat2,topo2),shading interp, latratio, caxis([-5 2])
% title('One-sixth degree interpolation')
% 
% region=[350 15 54 65];
% [topo1,lat1,lon1]=readtopo(region);
% [topo2,lat2,lon2]=readtopo(region,1/6);
% 
% subplot(2,2,3)
% pcolor(lon1,lat1,topo1),shading interp, latratio, caxis([-5 2])
% title('Original resolution Smith and Sandwell')
% subplot(2,2,4)
% pcolor(lon2,lat2,topo2),shading interp, latratio, caxis([-5 2])
% title('One-sixth degree interpolation')

figure,region=[-62 4 47 68];
[topo,lat,lon]=readtopo(region);
jpcolor(lon,lat,topo),latratio(58),hold on,
caxis([-5 3]),contour(lon,lat,topo,[0 0],'k')
xlabel('Longitude'),ylabel('Latitude')
title('Northern North Atlantic subset of Smith and Sandwell Topography')

if 0
    currentdir=pwd;
    cd([whichdir('jlab_license') '/figures'])
    print -dpng readtopo
    crop readtopo.png
    cd(currentdir)
end


function[]=readtopo_test

if isempty(which('jdata'))
     return
end

region=[-70 -35 50 70];
[topo1,lat1,lon1]=readtopo(region);
[long,latg]=meshgrid(lon1,lat1);
bool=inregion(region,latg,long);
reporttest('READTOPO strictly within region',allall(bool==1))

% if exist('extract_1m')==2  
%     [topo2,lat2,lon2]=extract_1m([region(3:4) region(1:2)]);
% 
%     [long,latg]=meshgrid(lon2,lat2);
%     bool=inregion(region,latg,long);
%     topo2=topo2(bool);
%     topo2=reshape(topo2,2492,2100);
% 
%     reporttest('READTOPO matches EXTRACT_1M',aresame(topo1,topo2/1000))
% end

%[topo1,lat1,lon1]=readtopo([-180 180 -80.738 -80.738]);

[topo,lat,lon,offset]=readtopo([0 0 80.738 80.738]);
reporttest('READTOPO first point sanity check',length(topo)==1&&offset==0)

%Note extract_1m gives 43200, so have have missed the northernmost stripe
%[topo,lat,lon]=extract_1m([80.738 80.738 0 0 ]);

[topo,lat,lon,offset]=readtopo([360-1/60 360-1/60 -80.738 -80.738]);
reporttest('READTOPO last point sanity check',length(topo)==1&&offset==(746496000-2))

%Note extract_1m gives 43200, so have have missed the northernmost stripe
%[topo,lat,lon]=extract_1m([-80.738 -80.738 0 0 ]);

region=[-65 -40 54 65];
[topo1,lat1,lon1,offset]=readtopo(region);
region=[deg360_internal(-65) deg360_internal(-40) 54 65];
[topo2,lat2,lon2,offset]=readtopo(region);
reporttest('READTOPO 180 and 360 degree formats match',aresame(topo1,topo2))
%figure,pcolor(lon1,lat1,topo1),shading flat

region=[165 -170 54 65];
[topo1,lat1,lon1,offset]=readtopo(region);
region=[deg360_internal(165) deg360_internal(-170) 54 65];
[topo2,lat2,lon2,offset]=readtopo(region);
reporttest('READTOPO 180 and 360 degree formats match, overlapping dateline',aresame(topo1,topo2))
%figure,pcolor(lon1,lat1,topo1),shading flat

region=[350 15 54 65];
[topo1,lat1,lon1,offset]=readtopo(region);
region=[deg180_internal(350) deg180_internal(15) 54 65];
[topo2,lat2,lon2,offset]=readtopo(region);
reporttest('READTOPO 180 and 360 degree formats match, overlapping prime meridian',aresame(topo1,topo2))
%figure,pcolor(lon1,lat1,topo1),shading flat

[topo,lat,lon,offset]=readtopo([0 360 80.738 80.738]);
bool(1)=length(topo)==21600;
bool(2)=offset==0;
bool(3)=lon(1)==1/120;
bool(4)=lon(end)==360-1/120;
bool(5)=aresame(80.738-lat,(1/120)/(1/cosd(80.738)),2e-5);
reporttest('READTOPO northernmost stripe test',allall(bool))

[topo,lat,lon,offset]=readtopo([0 360 -80.738 -80.738]);
bool(1)=length(topo)==21600;
bool(2)=offset+length(topo)*2==746496000;
bool(3)=lon(1)==1/120;
bool(4)=lon(end)==360-1/120;
bool(5)=aresame(-80.738-lat,-(1/120)/(1/cosd(80.738)),2e-5);
reporttest('READTOPO southernmost stripe test',allall(bool))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Helper functions below here, to be self-contained

function[varargout]=deg180_internal(varargin)
%DEG180  Converts degrees to the range [-180,180].
%
%   [TH1,TH2,...THN]=DEG180(TH1,TH2,...THN) converts the input
%   angles, which are measured in degrees, to the range [-180, 180].
%
%   See also DEG360.
%
%   'deg180 --t' runs a test.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007-2009 J.M. Lilly --- type 'help jlab_license' for details
 
varargout=varargin;
for i=1:nargin;
    varargout{i}=deg360_internal(varargin{i});
    bool=varargout{i}>180;
    varargout{i}(bool)=varargout{i}(bool)-360;
end

function[varargout]=deg360_internal(varargin)
%DEG360  Converts degrees to the range [0, 360].
%
%   [TH1,TH2,...THN]=DEG360(TH1,TH2,...THN) converts the input
%   angles, which are measured in degrees, to the range [0, 360].
%
%   More precisely, the output range is defined to include 0, but exclude
%   360, so DEG360(360) is defined as 360-EPS where EPS is 1e-10.
%
%   See also DEG180.
%
%   'deg360 --t' runs a test.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2013 J.M. Lilly --- type 'help jlab_license' for details
 
varargout=varargin;
for i=1:nargin;
    temp=varargin{i};
    temp(temp==360)=(360-1e-10); %Correct for 360 
    varargout{i}=mod(temp,360);

    bool=~isfinite(varargin{i});
    varargout{i}(bool)=varargin{i}(bool);
end


