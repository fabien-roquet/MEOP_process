function[mat,xmid,ymid,index]=twodhist(varargin)
%TWODHIST  Two-dimensional histogram.
%   __________________________________________________________________
%
%   *|* drifters_hist.png --- Histogram of the global drifter dataset.  
%   Type 'jhelp twodhist' to view this image. *|*
%   __________________________________________________________________
%
%   MAT=TWODHIST(X,Y,XBIN,YBIN) where X and Y are arrays of the same
%   length, creates a two-dimensional histogram MAT with bin edges
%   specified by XBIN and YBIN. 
%
%   If XBIN and YBIN are length N and M, respectively, then MAT is of
%   size M-1 x N-1.  XBIN and YBIN must be monotonically increasing. 
%
%   [MAT,XMID,YMID]=TWODHIST(...) optionally returns the midpoints XMID
%   and YMID of the bins.
%
%   TWODHIST, TWODSTATS, and TWODMED are three related functions for 
%   computing statistics as a function two variables using very fast
%   algorithms that avoid any loops through efficient use of indexing.
%
%   X and Y can also be cell arrays of numerical arrays, in which case 
%   all data values are concatented prior to finding the histogram.
%   __________________________________________________________________
%
%   Automatic bin calculation
%
%   TWODHIST can compute appropriate bins internall.
%
%   [MAT,XMID,YMID]=TWODIST(X,Y,N) uses N bins in the X and Y directions,
%   linearly spaced between the minimum and maximum values, and returns the
%   bin midpoints in XMID and YMID.  MAT is N-1 x N-1.
%
%   [MAT,XMID,YMID]=TWODIST(X,Y,[XMIN XMAX],[YMIN YMAX],N) similarly uses N
%   bins, linearly spaced between the designated X and Y values.  
%   __________________________________________________________________
%
%   Additional output
%   
%   [MAT,XMID,YMID,INDEX]=TWODHIST(X,Y,...) also returns INDEX, an array of
%   the same size as X and Y giving the index into the matrix MAT
%   corresponding to each (X,Y) data point. 
%   __________________________________________________________________
% 
%   Algorithms
%
%   TWODHIST uses a very fast and exact algorithm which is particularly 
%   efficient for large arrays.
%  
%   The (X,Y) data points are sorted according to their X-values, which 
%   determine a column within MAT, and their Y-values, which determine 
%   a row within MAT.  The historgram is formed looplessly, using indexing.
%
%   TWODHIST(...,'slow') uses a slow algorithm which uses less memory.  By
%   default, TWODHIST uses the fast but memory-intensive algorithm.  Use
%   this flag if you get an out-of-memory error.  
%   __________________________________________________________________
%
%   See also TWODMED, TWODSTATS.
%
%   'twodhist --f' generates the sample figure shown above.
%   'twodhist --t' runs a test.
%
%   Usage: [mat,xmid,ymid]=twodhist(x,y,N);
%          [mat,xmid,ymid]=twodhist(x,y,[xmin xmax],[ymin ymax],N);
%          mat=twodhist(x,y,xbin,ybin);
%          [mat,xmid,ymid]=twodhist(x,y,xbin,ybin);
%          [mat,xmid,ymid,index]=twodhist(x,y,xbin,ybin);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2015 J.M. Lilly --- type 'help jlab_license' for details    

  
if strcmpi(varargin,'--t')
   twodhist_test;return
elseif strcmpi(varargin,'--f')
   twodhist_figure;return
end

str='fast';
xdata=varargin{1};
ydata=varargin{2};

if iscell(xdata)
    [xdata,ydata]=cell2col(xdata,ydata);
end
if ~isreal(xdata)||~isreal(ydata)
    error('X and Y must be real-valued.');
end
if ~aresame(size(xdata),size(ydata))
     error('X and Y should have the same size.')
end
vcolon(xdata,ydata);

bool=isfinite(xdata)&isfinite(ydata);
xdata=xdata(bool);
ydata=ydata(bool);

if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
end

if length(varargin)==4
    xbin=varargin{3};
    ybin=varargin{4};
elseif length(varargin)==3
    N=varargin{3};
    xbin=linspace(minmin(xdata),maxmax(xdata),N);
    ybin=linspace(minmin(ydata),maxmax(ydata),N);
elseif length(varargin)==5
    N=varargin{5};
    xbin=linspace(varargin{3}(1),varargin{3}(2),N);
    ybin=linspace(varargin{4}(1),varargin{4}(2),N);
end
    
xbin=xbin(:);
ybin=ybin(:);
if any(diff(xbin)<0)
  error('XBIN must be monotonically increasing')
end
if any(diff(ybin)<0)
  error('YBIN must be monotonically increasing')
end

%Exclude points which are obviously outside of the domain
bool=xdata<xbin(end)&xdata>xbin(1)&ydata<ybin(end)&ydata>ybin(1);
xdata=xdata(bool);
ydata=ydata(bool);

if ~isempty(xdata)
    if ~isempty(strfind(str,'fast'))
      [mat,index]=twodhist_fast(xdata,ydata,xbin,ybin);
    else
      mat=twodhist_slow(xdata,ydata,xbin,ybin);
    end
else
    disp('Warning: No valid data in specified region.')
    mat=0*oprod(ybin(1:end-1),xbin(1:end-1)); 
end

if nargout>1
  xmid=(xbin+vshift(xbin,1,1))./2;
  xmid=xmid(1:end-1);
end
if nargout>2
  ymid=(ybin+vshift(ybin,1,1))./2;
  ymid=ymid(1:end-1);
end

function[mat,index]=twodhist_fast(xdata,ydata,xbin,ybin)

[xnum,xi,xmid]=bindata(xbin,xdata);
[ynum,yi,ymid]=bindata(ybin,ydata);

mat=zeros([length(ybin)-1,length(xbin)-1]);
index=nan*zeros(size(xdata));

nani=(~isnan(xnum)&~isnan(ynum));

if sum(nani(:))>0
    index(nani)=sub2ind([length(ybin)-1,length(xbin)-1],ynum(nani),xnum(nani));
    [indexsorted,sorter]=sort(index(nani));
    
    [L,ia]=blocklen(indexsorted);
    mat(indexsorted(ia))=L(ia);
end





function[mat]=twodhist_slow(xdata,ydata,xbin,ybin)
mat=0*oprod(ybin,xbin); 
[xbinb,ybinb]=vshift(xbin,ybin,1,1);
for i=1:length(xbin)
   for j=1:length(ybin)
         mat(j,i)=length(find(xdata>xbin(i)&xdata<=xbinb(i)&...
			      ydata>ybin(j)&ydata<=ybinb(j)));
   end
end
mat=mat(1:end-1,:);
mat=mat(:,1:end-1); 

function[]=twodhist_test
L=10000;
xdata=3*abs(rand(L,1));
ydata=3*abs(rand(L,1));
xbin=(0:.1:2);
ybin=(0:.2:2);
tic;
mat1=twodhist(xdata,ydata,xbin,ybin,'fast');
dt1=toc;
tic
mat2=twodhist(xdata,ydata,xbin,ybin,'slow');
dt2=toc;
bool=aresame(mat1,mat2,1e-10);
reporttest('TWODHIST fast vs. slow algorithm',bool)
disp(['TWODHIST hist algorithm was ' num2str(dt2./dt1) ' times faster than direct algorithm.'])

xdata=-3*abs(rand(L,1));
ydata=-3*abs(rand(L,1));
xbin=(-2:.1:0);
ybin=(-2:.2:0);
mat1=twodhist(xdata,ydata,xbin,ybin,'fast');
mat2=twodhist(xdata,ydata,xbin,ybin,'slow');
bool=aresame(mat1,mat2,1e-10);
reporttest('TWODHIST fast vs. slow algorithm, negative bins',bool)

xdata=randn(L,1);
ydata=randn(L,1);
xbin=(-2:.1:2);
ybin=(-2:.2:2);
mat1=twodhist(xdata,ydata,xbin,ybin,'fast');
mat2=twodhist(xdata,ydata,xbin,ybin,'slow');
bool=aresame(mat1,mat2,1e-10);
reporttest('TWODHIST fast vs. slow algorithm, crossing zero',bool)




function[]=twodhist_figure

if isempty(which('drifters.mat'))
    disp('Sorry, TWODHIST can''t find DRIFTERS.MAT.')
    return
end

disp('This make take a few minutes...')

load drifters,use drifters
[mat,xmid,ymid]=twodhist(lon,lat,-180.5:180.5,-89.5:89.5);

figure,jpcolor(xmid,ymid,log10(mat))
axis([-180 180 -70 80]),latratio(30),topoplot
xlabel('Longitude'),ylabel('Latitude')
title('Histogram of Data from the Global Surface Drifter Dataset')
caxis([1.5 3.5]),h=colorbar('EastOutside');
axes(h),ylabel('Log10 Number of Observations')

if 0
    currentdir=pwd;
    cd([whichdir('jlab_license') '/figures'])
    print -dpng drifters_hist
    crop drifters_hist.png
    cd(currentdir)
end

