function[h,indexout]=ellipseplot(varargin)
%ELLIPSEPLOT  Plot ellipses.
%
%   ELLIPSEPLOT(K,L,TH,Z) plots an ellipse of amplitude K, linearity L, and
%   orientation TH at complex-valued location Z=X+iY.
%
%   K, L, and TH are shorthand for the symbols KAPPA, LAMBDA, and THETA.
%   K and L are related to the semi-axes A and B by K=SQRT((A^2+B^2)/2)
%   and L=(A^2-B^2)/(A^2+B^2).  See Lilly and Gascard (2006) for details. 
%
%   Multiple ellipses are plotted if K, L, TH, and Z, are arrays of the
%   same size.
%
%   ELLIPSEPLOT draws each column of the input arrays with a different 
%   color, cycling through the default line colors.
%
%   H=ELLIPSEPLOT(...) outputs an array of line handles.
%
%   ELLIPSEPLOT calls ELLCURVES to actually compute the elliptical curves. 
%   ______________________________________________________________________
%
%   Aspect ratio
%
%   ELLIPSEPLOT(...,AR) with AR=[XAR YAR] multiplies the X-signal and
%   Y-signal by XAR and YAR, respectively, for plotting purposes. The
%   aspect ratio is set such that circles appear circular.
%
%   It's okay for AR to have three entries, AR=[XAR YAR ZAR], as is
%   output by "get(gca,'dataaspectratio')".  The last entry is ignored.
%
%   AR is optional and defaults to [1 1]. 
%
%   The following trailing options can occur in any order, as long as 
%   they are after K, L, TH, Z, and (if present) AR.
%   ______________________________________________________________________
%
%   Skipping ellipses
%
%   Frequently one does not wish to plot all the ellipses.  There are two
%   ways to accomplish this, as documented in this and the next section.
%
%   K, L, TH, and Z are column vectors, plot every SKIP-th ellipse with
%
%         ELLIPSEPLOT(K,L,TH,Z, ... ,'skip',SKIP).
%
%   This plots ellipses at indicies [SKIP:SKIP:LENGTH(K)-SKIP].
%
%   More generally, if the input field are matrices of the same size, each
%   having N dimensions, then SKIP can be an array with N elements.  Then
%   SKIP(1) indicates the SKIP number for the first matrix dimension, etc.
%   ______________________________________________________________________
%
%   Indexing ellipses
%   
%   To have more precise control over which ellipses are plotted, use
%   
%          ELLIPSEPLOT(K,L,TH,Z, ... ,'index',INDEX) 
%
%   which only plots the ellipses at the indicies indicated by INDEX. If
%   INDEX is empty, nothing happens.  
%
%   When K, L, TH, and Z are column vectors, INDEX is an array.  
% 
%   More generally, INDEX can be a cell array of N arrays, one for each 
%   dimension of the input matrices.
%
%   See PERIODINDEX for generating an index that skips every N periods.
%   ______________________________________________________________________
%
%   Specifying linestyles
%
%   ELLIPSEPLOT(K,L,TH,Z, ... ,STY) where STY is a single style string in 
%   LINESTYLE format, or a cell array of such strings, determines the line
%   styles for the ellipses.
%
%   With cell array input, the entries of STY is applied to the elements in
%   the cell arrays K{1}, K{2}, etc.  Otherwise, the entries in STY are 
%   applied to the multiple columns (if any) of the arrays K, L, TH and Z.
%
%   By default, the default line styles are used, so corresponding to 
%   STY{1}='b'; STY{2}='g'; STY{3}='r'; and so forth.
%
%   For example, STY{1}='2b--'; STY{2}='r-.'; will alternate between thick
%   blue dashed lines and thin red dash-dotted lines.
%
%   See LINESTYLE for more details. 
%   ______________________________________________________________________
% 
%   Additional options
%
%   ELLIPSEPLOT(K,L,TH,Z, ... ,'m_map') will work with Rich Pawlowicz's 
%   M_MAP package. In this case Z should be of the form LON+SQRT(-1)*LAT,
%   or a cell array of numeric arrays having this form.
%
%   ELLIPSEPLOT(K,L,TH,Z, ... ,'phase',PHI) optionally draws a small line, 
%   like a clock hand, to indicate the ellipse phase PHI.
%
%   ELLIPSEPLOT(K,L,TH,Z, ... ,'npoints',N) plots ellipses with N points
%   along the circumference.  The default value is 32.  Use N=16 or 64 
%   for faster plotting or for smoother, more well-defined ellipses.
%
%   ELLIPSEPLOT(K,L,TH,Z, ... ,'scale',S) sets the ellipse scale for a 
%   plot in which Z is of the form LON+SQRT(-1)*LAT.  Here K should be
%   in units kilometers, and the resulting ellipses will be plotted at S
%   times their actual size, with S=1 corresponding to actual size.
%
%   [H,INDEX]=ELLIPSEPLOT(K,L,TH,Z,...) also returns the index INDEX into
%   the ellipse locations that were actually mapped.  This is particularly 
%   useful with the 'skip' input argument, and allows the handles H to 
%   be manipulated later, for example, to set color based on K values.
%   ______________________________________________________________________
%
%   Cell array input 
%
%   ELLIPSEPLOT also works if the input arguments are cell arrays. 
%
%   Specifically, if K, L, TH, and Z are all cell arrays of N different 
%   numerical arrays, then ELLIPSEPLOT will act on each of them in turn.
%
%   In this case SIZE(K{M}), SIZE(L{M}), SIZE(TH{M}), and SIZE(Z{M}) must 
%   all be identical for M=1,2,... N.
%   ______________________________________________________________________
%
%   See also PERIODINDEX, LINECOLOR, ELLCURVES.
%
%   'ellipseplot --f' makes some test figures.
%
%   Usage: ellipseplot(k,l,th,z)
%          ellipseplot(k,l,th,z,ar)
%          ellipseplot(k,l,th,z,'axis')
%          ellipseplot(k,l,th,z,ar,'phase',phi)
%          ellipseplot(k,l,th,z,ar,'npoints',64)
%          ellipseplot(k,l,th,z,ar,'index',index)
%          ellipseplot(k,l,th,z,ar,'skip',5)
%          ellipseplot(k,l,th,z,'2r--')
%   ______________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2014 J.M. Lilly --- type 'help jlab_license' for details        

if strcmpi(varargin{1},'--f')
  ellipseplot_fig;return
end

%/********************************************************
%Sort out input arguments
na=length(varargin);

k=varargin{1};
l=varargin{2};
th=varargin{3};
sty{1}='b';
sty{2}='g';
sty{3}='r';
sty{4}='c';
sty{5}='m';
sty{6}='y';
sty{7}='k';

str='matlab';

naold=na+1;
index=[];
skip=[];
phi=[];
x=[];
scale=[];
npoints=32;

while naold~=na
    naold=na;   
    if ischar(varargin{end})
        if length(varargin{end})<3
            clear sty
            sty{1}=varargin{end};
        elseif ~strcmpi(varargin{end}(1:3),'axi')&&~strcmpi(varargin{end}(1:3),'mat')&&...
                ~strcmpi(varargin{end}(1:3),'m_m')
            clear sty
            sty{1}=varargin{end};
        else
            str=varargin{end};
        end       
        varargin=varargin(1:end-1);
        na=na-1;     
    end
       
    if iscell(varargin{end})&&~strcmpi(varargin{end-1}(1:3),'ind')
       sty=varargin{end};
       varargin=varargin(1:end-1);
       na=na-1;  
    end
    
    if ischar(varargin{end-1})
      if strcmpi(varargin{end-1},'phase')
          phi=varargin{end};
      elseif strcmpi(varargin{end-1},'npoints')
          npoints=varargin{end};
      elseif strcmpi(varargin{end-1},'index')
          index=varargin{end};
      elseif strcmpi(varargin{end-1},'skip')
          skip=varargin{end};
      elseif strcmpi(varargin{end-1},'scale')
          scale=varargin{end};
      end
      na=na-2;
      varargin=varargin(1:end-2);
    end
end

if na>3
  x=varargin{4};
end

if na>4
  ar=varargin{5};
  if size(ar,2)~=2&&size(ar,2)~=3
    error('The aspect ratio AR must have length equal to 2 or 3.')
  end
  if length(ar)>3,warning('I was expecting AR to be length 2 or 3.'),end
  ar1=ar(1);
  ar2=ar(2);
else
  ar1=1;
  ar2=1;
end
%\********************************************************
if isempty(phi)
    if iscell(k)
        for i=1:length(k)
             phi{i,1}=nan*zeros(size(k{i}));
        end
    else
        phi=nan*zeros(size(k));
    end 
end

if isempty(x)
    if iscell(k)
        for i=1:length(k)
            x{i}=zeros(size(k{i}));
        end
    else
        x=zeros(size(k));
    end 
end


bhold=ishold;
set(gca,'visible','off')
hold on
storestate=get(gcf,'BackingStore');
set(gcf,'BackingStore','off')

%vsize(k,l,th,phi,x,index,skip,ar1,ar2,npoints,scale)
%index

clear h indexout
if iscell(k)
    for i=1:length(k)
        %vsize(k{i},l{i},th{i},phi{i},x{i},ar1,ar2,baxis,npoints,str,index,skip,sty)
        if ~isempty(k{i})
            [h{i},indexout{i}]=ellipseplot_loop(k{i},l{i},th{i},phi{i},x{i},ar1,ar2,npoints,str,index{i},skip,scale);
        else
            h{i}=[];
            indexout{i}=[];
        end
    end
else
    %vsize(k,l,th,phi,x)
    %npoints,str, index,skip,scale
    [h,indexout]=ellipseplot_loop(k,l,th,phi,x,ar1,ar2,npoints,str,index,skip,scale);
end

set(gca,'dataaspectratio',[ar1 ar2 1])
if iscell(h)
    for i=1:length(h)
        %I can deal with this by forming a long string like '2b','2b', etc
        %hindex=res((id-i)./length(sty))==0;
        %if ~isempty(hindex)
        linestyle(h{i},sty{mod(i-1,length(sty))+1});
        %end
    end
else
    for i=1:size(h,2)
        sty{mod(i-1,length(sty))+1};
        linestyle(h(:,i),sty{mod(i-1,length(sty))+1});
    end
end

if verLessThan('matlab','8.4.0')
    if iscell(h)
        h2=cell2col(h);
    else
        h2=h;
    end
    set(h2(~isnan(h2)),'visible','on')
else
    set(h,'visible','on')
end
set(gcf,'BackingStore',storestate)
if ~bhold
    hold off
end

if ~strcmpi(str(1:3),'m_m')
    set(gca,'box','on')
    set(gca,'visible','on')
end

if nargout==0
    clear h
end

function[h,indexfull]=ellipseplot_loop(k,l,th,phi,x,ar1,ar2,npoints,str,index,skip,scale)

if ~iscell(index)&&~isempty(index)
    temp=index;
    clear index;
    index{1}=temp;
end

if ~isempty(skip)
    for i=1:length(skip)
        index{i}=skip(i):skip(i):size(k,i)-skip(i);
    end
end

h=[];
indexfull=[1:length(k(:))]';
indexfull=reshape(indexfull,size(k));
if ~isempty(index)
    for i=1:length(index)
        if ~isempty(index{i})
            vindex(k,l,th,phi,x,indexfull,index{i},i);
        else
            [k,l,th,phi,x,indexfull]=vempty;
        end
    end
end

clear index


if ~isempty(k)
    id=vrep(1:size(k,2),size(k,1),1);
    vcolon(k,l,th,phi,x,indexfull,id);
%    k,l,th,phi,x,indexfull,id
    vswap(k,0,nan);
    index=find(~isnan(k));
    
    if ~isempty(scale)
        k=frac(k,frac(2*pi,360)*radearth)*scale;
    end
    
    %Make sure signal will be complex-valued
    l(l==0)=1e-10;
    if ~isempty(index)
        vindex(k,l,th,phi,x,id,indexfull,index,1);
        h=ellipseplot1(k,l,th,phi,x,ar1,ar2,npoints,str);
    end
else
    h=[];
end

function[h]=ellipseplot1(kappa,lambda,theta,phi,x,ar1,ar2,npoints,str)
%kappa,lambda,theta,phi,x,ar1,ar2,npoints

%New algorithm as of May 2014, much faster
if ~allall(isnan(phi));
    z=ellcurves(kappa,lambda,theta,phi,x,'npoints',npoints,'aspect',[ar1,ar2]);
    z(end+1,:)=x;%To plot the phase clock
else
    z=ellcurves(kappa,lambda,theta,x,'npoints',npoints,'aspect',[ar1,ar2]);
end

%figure,plot(z)
if strcmpi(str(1:3),'mat')
    h=plot(z,'visible','off');%set(gca,'dataaspectratio',[ar(:)' 1]);
elseif strcmpi(str(1:3),'m_m')
    h=m_plot(real(z),imag(z),'visible','off');
end


function[]=ellipseplot_fig

a=ones(10,1);
b=(1:10)'./10;
[k,l]=ab2kl(a,b);
th=linspace(0,pi,10)';
x=linspace(0,1,10)'*20;

figure,
ellipseplot(k,l,th,x,'phase',0*th,'npoints',16,'g')
title('Counterclockwise rotating ellipse becoming circle')
set(gca,'dataaspectratio',[1 1 1])



h=ellipseplot(k,l,th,x,'npoints',16,'g','ellipses');
lato=5:10:85;
lono=0;

phi=[0:1:1000]'/10;
kappa=100+0*phi;
lambda=phi./maxmax(phi);
theta=(pi/2).*phi./maxmax(phi);
z=ellsig(kappa,lambda,theta,phi);
figure
for i=1:9
    subplot(3,3,i)
    [lat,lon]=xy2latlon(z,lato(i),lono);
    plot(lon,lat,'r'),latratio(lato(i)),hold on,
    ellipseplot(kappa,lambda,theta,lono+sqrt(-1)*lato(i)+0*kappa,latratio(lato(i)),'skip',100,'scale',1);
    axis tight
    if i==2
        title('Checking plotted ellipse size and shape at various latitudes')
    end
    xlabel('Longitude'),ylabel('Latitude')
    %text(-1,lato(i)-1,['Latitude ' int2str(lato(i))])
end


% for i=1:9
%     figure
%     [lat,lon]=xy2latlon(z,lato(i),lono);
%     %m_proj('albers equal-area conic','lon',[min(lon) max(lon)],'lat',[min(lat) max(lat)]);
%     m_proj('miller','lon',[min(lon) max(lon)],'lat',[min(lat) max(lat)]);
%     m_grid;m_coast('line','color','k'),hold on  
%     m_plot(lon,lat,'r'),hold on,
%     ellipseplot(kappa,lambda,theta,lono+sqrt(-1)*lato(i)+0*kappa,latratio(lato(i)),'skip',100,'scale',1,'m_map');
% end

