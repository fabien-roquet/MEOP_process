function[mat,xmid,ymid,num,std]=twodstats(varargin)
%TWODSTATS  Mean, variance, and covariance of functions of two variables.
%   __________________________________________________________________
%
%   *|* drifters_mean.png --- Mean speed from the global drifter dataset.  
%   Type 'jhelp twodstats' to view this image. *|*
%   __________________________________________________________________
%
%   TWODSTATS computes the first- and second-order statistics of a function
%   of two variables in prescribed bins.  This function may either be a 
%   scalar-valued or vector-valued quantity at each point. 
% 
%   An example of a scalar-valued dataset is temperature as a function of
%   latitude and longitude. An example of a vector-valued dataset is wind
%   or current velocity as a function of latitude and longitude.
%
%   TWODSTATS, TWODHIST, and TWODMED are three related functions for 
%   computing statistics as a function two variables using very fast
%   algorithms that avoid any loops through efficient use of indexing. 
%   __________________________________________________________________
%  
%   Mean and standard deviation of a scalar-valued function 
%
%   MZ=TWODSTATS(X,Y,Z,XBIN,YBIN) where X, Y and Z are arrays of the same
%   length, forms the mean of Z over the XY plane.  
%  
%   If XBIN and YBIN are length N and M, respectively, then MZ is of 
%   size M-1 x N-1.  Bins with no data are assigned a value of NAN.
%
%   XBIN and YBIN must be monotonically increasing. 
%
%   MZ=TWODSTATS(X,Y,Z,N) uses N bins in the X and Y directions, linearly
%   spaced between the minimum and maximum values.  MZ is N-1 x N-1.
%
%   MZ=TWODSTATS(X,Y,Z,[XMIN XMAX],[YMIN YMAX],N) uses N bins, linearly
%   spaced between the designated X and Y values.  MZ is N-1 x N-1. 
%
%   X, Y, and Z can also be cell arrays of numerical arrays, in which case 
%   all data values are concatented prior to finding the statistics.
%   __________________________________________________________________
%
%   Additional output
%   
%   [MZ,XMID,YMID]=TWODSTATS(...) optionally returns the midpoints XMID
%   and YMID of the bins.
%
%   [MZ,XMID,YMID,NUMZ]=TWODSTATS(...) also returns the number of good
%   data points in each of the (X,Y) bins.  NUMZ is the same size as MZ.
%
%   [MZ,XMID,YMID,NUMZ,STDZ]=TWODSTATS(...) also returns the standard 
%   deviation of Z in the (X,Y) bins.  STDZ is the same size as MZ.
%   __________________________________________________________________
%  
%   Mean and covariance of a vector-valued function 
%   
%   TWODSTATS can also be used to analyze a function which contains more
%   than one value at each (X,Y) point.  
%
%   If Z represents a vector with K components, then Z should have the same
%   size as X and Y in all but its last dimension, which will be length K.
%
%   MZ=TWODSTATS(X,Y,Z,XBIN,YBIN) then returns MZ, containing the mean 
%   values of each component of Z in each bin.  If M and N are the lengths 
%   of XBIN and YBIN, MZ is of size M-1 x N-1 x K. 
%
%   [MZ,XMID,YMID,NUMZ,COVZ]=TWODSTATS(...) returns the full covariance
%   matrix COVZ in each of the bins.  As the covariance of Z is K x K, the 
%   size of the output matrix COVZ is M-1 x N-1 x K x K.
%   __________________________________________________________________
%
%   Algorithms
%
%   TWODSTATS uses a fast (exact) algorithm which is particularly efficient 
%   for large arrays.  TWODSTATS can be two orders of magnitude faster than
%   the obvious way of sorting the data in bins using explicit loops.   
%
%   The values of Z are sorted into bins according to the (X,Y) values, 
%   with bin edges specified by XBIN and YBIN. The mean of all finite 
%   values of Z in each bin is then computed looplessly, using indexing.
%
%   TWODSTATS(...,'slow') uses a slow algorithm which uses less memory.  By
%   default, TWODSTATS uses a fast but memory-intensive algorithm.  Use 
%   this flag if you get an out-of-memory error.  
%   __________________________________________________________________
%   
%   See also TWODHIST, TWODMED.
%
%   'twodstats --t' runs a test.
%   'twodstats --f' generates the sample figure shown above.
%
%   Usage: mz=twodstats(x,y,z,N);
%          mz=twodstats(x,y,z,[xmin xmax],[ymin ymax],N);
%          mz=twodstats(x,y,z,xbin,ybin);
%          [mz,xmid,ymid]=twodstats(x,y,z,xbin,ybin);
%          [mz,xmid,ymid,numz]=twodstats(x,y,z,xbin,ybin);
%          [mz,xmid,ymid,numz,stdz]=twodstats(x,y,z,xbin,ybin);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2015 J.M. Lilly --- type 'help jlab_license' for details    

%   You can use TWODSTATS for fast binning of data over the plane.  For
%   the case in which Z is so sparsely distributed over X and Y, such 
%   that no bins will have more than one entry, the mean in each bin
%   is just the value of the data point in the bin.  

if ischar(varargin{1})
    if strcmpi(varargin{1},'--t')
        twodstats_test;return
    elseif strcmpi(varargin{1},'--t')
        twodstats_fig;return
    end
end

str='fast';
xdata=varargin{1};
ydata=varargin{2};
zdata=varargin{3};

if iscell(xdata)
    [xdata,ydata,zdata]=cell2col(xdata,ydata,zdata);
end
if ~aresame(size(xdata),size(ydata))
     error('X and Y should have the same size.')
end
if ~isreal(xdata)||~isreal(ydata)||~isreal(zdata)
    error('X, Y, and Z must be real-valued.');
end

vcolon(xdata,ydata);
K=length(zdata(:))/length(xdata);
zdata=reshape(zdata,length(xdata),K);
  

bool=isfinite(xdata)&isfinite(ydata)&isfinite(sum(zdata,2));
xdata=xdata(bool);
ydata=ydata(bool);
zdata=zdata(bool,:);

if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
end

if length(varargin)==5
    xbin=varargin{4};
    ybin=varargin{5};
elseif length(varargin)==4
    N=varargin{4};
    xbin=linspace(minmin(xdata),maxmax(xdata),N);
    ybin=linspace(minmin(ydata),maxmax(ydata),N);
elseif length(varargin)==6
    N=varargin{6};
    xbin=linspace(varargin{4}(1),varargin{4}(2),N);
    ybin=linspace(varargin{5}(1),varargin{5}(2),N);
end
    
    
if nargin==5
  str='fast';
end

vcolon(xbin,ybin);
if any(diff(xbin)<0)
  error('XBIN must be monotonically increasing')
end
if any(diff(ybin)<0)
  error('YBIN must be monotonically increasing')
end
    
if nargout>4
    stdflag=1;
else
    stdflag=0;
end

%Exclude points which are obviously outside of the domain
bool=xdata<xbin(end)&xdata>xbin(1)&ydata<ybin(end)&ydata>ybin(1);
xdata=xdata(bool);
ydata=ydata(bool);
zdata=zdata(bool,:);

if ~isempty(zdata)
    if ~isempty(strfind(str,'fast'))
      [mat,num,std]=twodstats_fast(xdata,ydata,zdata,xbin,ybin,stdflag);
    elseif ~isempty(strfind(str,'slow'))
      [mat,num,std]=twodstats_slow(xdata,ydata,zdata,xbin,ybin,stdflag);
    end
else
    disp('Warning: Data contains only NANs and / or INFs.')
    mat=0*oprod(ybin(1:end-1),xbin(1:end-1)); 
    num=mat;
    std=mat;
end


if nargout>1
  xmid=(xbin+vshift(xbin,1,1))./2;
  xmid=xmid(1:end-1);
end
if nargout>2
  ymid=(ybin+vshift(ybin,1,1))./2;
  ymid=ymid(1:end-1);
end

index=find(num==0);
if ~isempty(index)
    mat(index)=nan;
    if stdflag
        %for i=1:size(std,3)...
        std(index)=nan;
    end
end

function[mat,num,covmat]=twodstats_fast(xdata,ydata,zdata,xbin,ybin,stdflag)
%vsize(xdata,ydata,zdata,xbin,ybin,stdflag)
[num,mat,covmat]=vzeros((length(ybin)-1)*(length(xbin)-1),1,'nan');
mat=vrep(mat,size(zdata,2),2);
covmat=vrep(vrep(covmat,size(zdata,2),2),size(zdata,2),3);
xnum=bindata(xbin,xdata);
ynum=bindata(ybin,ydata);

index=sub2ind([length(ybin)-1,length(xbin)-1],ynum,xnum);

%figure,plot(xnum,ynum,'.')
%vsize(index,xdata,ydata,zdata)


if ~isempty(index)
    [index,sorter]=sort(index);
    zdata=double(zdata(sorter,:));
    [L,ia,ib,numblock]=blocklen(index);
    
    num(index(ia))=L(ia);
    vswap(num,0,nan);
     
    if 1
        %Cumsum has problem: Force double precision for good results
        cumsumzdata=cumsum(zdata,1);
        mat(index(ia),:)=cumsumzdata(ib,:)-cumsumzdata(ia,:)+zdata(ia,:);
    else
        %For testing purposes
        for i=1:length(ia)
            mat(index(ia(i)),:)=sum(zdata(ia(i):ib(i)),:);
        end
    end
    mat=mat./vrep(num,size(mat,2),2);

    if stdflag
        %This trickery simply puts the mean back where it belongs in the column
        blockmean=zeros(length(index),size(zdata,2)); 
        blockmean(1,:)=mat(index(ia(1)),:);
        blockmean(ia(2:end),:)=mat(index(ia(2:end)),:)-mat(index(ia(1:end-1)),:);
        blockmean=cumsum(blockmean,1);
        
        zprime=zdata-blockmean;
        if size(zdata,2)==1        
            zvar=abs(zprime).^2;
        else
            zprime2=conj(permute(zprime,[1 3 2]));
            zprime=vrep(zprime,size(zdata,2),3);
            zprime2=vrep(zprime2,size(zdata,2),2);
            zvar=zprime.*zprime2;
        end
        
        if 1
            %Cumsum has problem: Force double precision for good results
            cumsumzvar=cumsum(double(zvar),1);
            covmat(index(ia),:,:)=cumsumzvar(ib,:,:)-cumsumzvar(ia,:,:)+zvar(ia,:,:);
        else
            %Just so you know what the above line actually does
            for i=1:length(ia)
                covmat(index(ia(i)),:,:)=sum(zvar(ia(i):ib(i)),:,:);
            end
        end
        
        covmat=covmat./vrep(vrep(num,size(covmat,2),2),size(covmat,3),3);
                
        if size(zdata,2)==1 
            covmat=sqrt(covmat);
        end
            
    end
end
vswap(num,nan,0);
mat=reshape(mat,length(ybin)-1,length(xbin)-1,size(zdata,2));
num=reshape(num,length(ybin)-1,length(xbin)-1);
covmat=reshape(covmat,length(ybin)-1,length(xbin)-1,size(zdata,2),size(zdata,2));


function[mat,num,std]=twodstats_slow(xdata,ydata,zdata,xbin,ybin,stdflag)
vcolon(xdata,ydata,zdata,xbin,ybin);
index=find(isfinite(xdata)&isfinite(ydata)&isfinite(zdata));
vindex(xdata,ydata,zdata,index,1);

mat=0*oprod(ybin,xbin); 
num=0*oprod(ybin,xbin); 
    
if stdflag 
    std=0*oprod(ybin,xbin); 
else
    std=[];
end

[xbinb,ybinb]=vshift(xbin,ybin,1,1);
for i=1:length(xbin)
   for j=1:length(ybin)
         index=find(xdata>xbin(i)&xdata<=xbinb(i)&ydata>ybin(j)&ydata<=ybinb(j));
         if ~isempty(index)
             mat(j,i)=vmean(zdata(index),1);
             num(j,i)=length(index);
             if stdflag
                 std(j,i)=vstd(zdata(index),1);
             end
         end
   end
end
mat=mat(1:end-1,:);
mat=mat(:,1:end-1); 
num=num(1:end-1,:);
num=num(:,1:end-1);
    
if stdflag
    std=std(1:end-1,:);
    std=std(:,1:end-1);
end

function[]=twodstats_test
twodstats_test1;
twodstats_test2;

function[]=twodstats_test1
load npg2006
use npg2006
lono=[-19:.05:-16.5];
lato=[42.8:.05:44.3];


tic;
[mat1,xmid,ymid,num1,std1]=twodstats(lon,lat,t,lono,lato,'fast');
dt1=toc;
tic
[mat2,xmid,ymid,num2,std2]=twodstats(lon,lat,t,lono,lato,'slow');
dt2=toc;

%index=isfinite(mat1)&isfinite(mat2);
index=1:length(mat1(:));
%size(find(index))

%figure,pcolor(mat1),shading flat
%figure,pcolor(mat2),shading flat

bool1=aresame(mat1(index),mat2(index),1e-10);
bool2=aresame(std1(index),std2(index),1e-10);
bool3=aresame(num1(index),num2(index),1e-10);
reporttest('TWODSTATS fast vs. slow algorithm for npg2006, mean',bool1)
reporttest('TWODSTATS fast vs. slow algorithm for npg2006, std',bool2)
reporttest('TWODSTATS fast vs. slow algorithm for npg2006, num',bool3)
disp(['TWODSTATS fast algorithm was ' num2str(dt2./dt1) ' times faster than direct algorithm.'])


function[]=twodstats_test2
L=10000;
xdata=3*abs(rand(L,1));
ydata=3*abs(rand(L,1));
zdata=randn(L,1);
xbin=(0:.1:2);
ybin=(0:.2:2);
tic;
[mat1,xmid,ymid,num1,std1]=twodstats(xdata,ydata,zdata,xbin,ybin,'fast');
dt1=toc;
tic
[mat2,xmid,ymid,num2,std2]=twodstats(xdata,ydata,zdata,xbin,ybin,'slow');
dt2=toc;
bool1=aresame(mat1,mat2,1e-10);
bool2=aresame(std1,std2,1e-10);
bool3=aresame(num1,num2,1e-10);
reporttest('TWODSTATS fast vs. slow algorithm, mean',bool1)
reporttest('TWODSTATS fast vs. slow algorithm, std',bool2)
reporttest('TWODSTATS fast vs. slow algorithm, num',bool3)
disp(['TWODSTATS fast algorithm was ' num2str(dt2./dt1) ' times faster than direct algorithm.'])

xdata=-4*rand(L,1);
ydata=-5*rand(L,1);
zdata(1:round(L/10):end)=nan;
xbin=(-2:.1:0);
ybin=(-3:.2:0);
tic;[mat1,xmid,ymid,num1,std1]=twodstats(xdata,ydata,zdata,xbin,ybin,'fast');
dt1=toc;
tic;[mat2,xmid,ymid,num2,std2]=twodstats(xdata,ydata,zdata,xbin,ybin,'slow');
dt2=toc;
bool1=aresame(mat1,mat2,1e-10);
bool2=aresame(std1,std2,1e-10);
bool3=aresame(num1,num2,1e-10);
reporttest('TWODSTATS fast vs. slow algorithm, with negative bins and NANs, mean',bool1)
reporttest('TWODSTATS fast vs. slow algorithm, with negative bins and NANs, std',bool2)
reporttest('TWODSTATS fast vs. slow algorithm, with negative bins and NANs, num',bool3)

disp(['TWODSTATS fast algorithm was ' num2str(dt2./dt1) ' times faster than direct algorithm.'])


function[]=twodstats_fig

if isempty(which('drifters.mat'))
    disp('Sorry, TWODSTATS can''t find DRIFTERS.MAT.')
    return
end

disp('This make take a few minutes...')

figure
load drifters,use drifters
[mat,xmid,ymid]=twodstats(lon,lat,cellabs(cv),-180.5:180.5,-89.5:89.5);
jpcolor(xmid,ymid,mat)
axis([-180 180 -70 80]),latratio(30),topoplot
xlabel('Longitude'),ylabel('Latitude')
title('Mean Speed from the Global Surface Drifter Dataset')
caxis([0 60]),h=colorbar('EastOutside');
axes(h),ylabel('Mean Speed (cm/s)')
echo off


if 0
    currentdir=pwd;
    cd([jlab_path '/figures'])
    print -dpng drifters_mean
    crop drifters_mean.png
    cd(currentdir)
end



