function[varargout]=cellchunk(varargin)
%CELLCHUNK  Converts cell array data into uniform length 'chunks'.
%
%   Y=CELLCHUNK(X,L), where X is a cell array of variable-length numeric
%   arrays, extracts non-overlapping 'chunks' of the uniform length L and
%   returns these in Y, a matrix having L rows.
%   
%   Successive chunks of length L data in X are become successive columns
%   in Y, regardless of the cell to which they initially belonged.  
%   Residuals of length less then L at the end of each cell are discarded.
%
%   For example:
%
%       x{1}=[1 2 3 4]'; x{2}=[5 6 7]'; x{3}=[8 9]'; 
%
%       cellchunk(x,2) =  [1 3 5 8]
%                         [2 4 6 9]
%
%   This is useful, for example, in examining spectra from uniform-length 
%   intervals in Lagragnian data such as FLOATS.MAT or DRIFTERS.MAT.
%
%   If there is a remainder R in how many times L goes into the length of
%   the time series, FLOOR(R/2) points will be thrown away from the
%   beginning of the time series and CEIL(R/2) points from the end. 
%
%   [Y1,Y2,...,YN]=CELLCHUNK(X1,X2,...,XN,L) also works, where the XN are
%   all cell array of the same size.
% 
%   CELLCHUNK with no output arguments overwrite the original named input
%   variables. 
%   __________________________________________________________________
%   
%   Overlap
%
%   CELLCHUNK(...,L,'overlap') instead outputs chunks of length L with a 
%   50% overlap.  That is, successive columns of the output overlap by L/2.
%
%   In this case, remainders are discarded from the ends of the time 
%   series, instead of being split between the beginning and the end as in
%   the non-overlapping case.  
%   __________________________________________________________________   
%
%   See also INERTIALCHUNK.
%
%   'cellchunk --t' runs a test.
%
%   Usage: y=cellchunk(x,L);
%          [y1,y2,y3]=cellchunk(x1,x2,x3,L);
%          [y1,y2,y3]=cellchunk(x1,x2,x3,L,'overlap');
%          cellchunk(x1,x2,x3,L);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(varargin{1}, '--t')
    cellchunk_test,return
end

if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
else
    str='nooverlap';
end

L=varargin{end};
varargin=varargin(1:end-1);
na=length(varargin);
for j=1:na
    if strcmpi(str(1:3),'ove')
        varargout{j}=cellchunk_overlap(L,varargin{j});
    else
        varargout{j}=cellchunk_nooverlap(L,varargin{j});
    end
end

eval(to_overwrite(na));


function[x]=cellchunk_nooverlap(L,x)
for i=1:length(x)
    N=length(x{i});
    M=mod(N,L);
    index=1+floor(M/2):N-ceil(M/2);
    x{i}=vindex(x{i},index,1);
end
%x=vcellcat(x);  %Replaced vcellcat with cell2col
x=nonnan(cell2col(x));


x=reshape(x,L,length(x)/L);


function[x]=cellchunk_overlap(L,x)

%Create an index number
num=x;
numo=0;
for i=1:length(x)
    num{i}=[1:length(num{i})]'+numo;
    numo=num{i}(end)-1;
end

%First half, length L starting at 1
x1=x;
num1=num;
for i=1:length(x)
    N=length(x{i});
    M=mod(N,L);
    index=1:N-M;
    [x1{i},num1{i}]=vindex(x{i},num{i},index,1);
end
[x1,num1]=cell2col(x1,num1);
x1=nonnan(x1);
num1=nonnan(x1);
%[x1,num1]=vcellcat(x1,num1);
x1=reshape(x1,L,length(x1)/L);
num1=reshape(num1,L,length(num1)/L);
num1=num1(1,:);

%Second half, length starting at L/2
x2=x;
num2=num;
for i=1:length(x)
    %Toss out floor(L/2)+1 dat points
    x{i}=x{i}(floor(L/2)+1:end);
    num{i}=num{i}(floor(L/2)+1:end);
    N=length(x{i});
    M=mod(N,L);
    index=1:N-M;
    [x2{i},num2{i}]=vindex(x{i},num{i},index,1);
end
[x2,num2]=cell2col(x2,num2);
x2=nonnan(x2);
num2=nonnan(x2);
%[x2,num2]=vcellcat(x2,num2);
x2=reshape(x2,L,length(x2)/L);
num2=reshape(num2,L,length(num2)/L);
num2=num2(1,:);

%What we're doing here is interleaving...
%but, we have to take care that the columns are all in the correct order,
%which is why we go through all this trouble with NUM.
x=[x1 x2];
num=[num1 num2];
%figure,plot(num)
[sorted,sorter]=sort(num);
%figure,plot(sorted)
x=x(:,sorter);





function[]=cellchunk_test
 
load ebasnfloats
use ebasnfloats

[numc,latc,lonc]=cellchunk(num(3),lat(3),lon(3),200);
dnum=numc(1,2:end)-numc(end,1:end-1);%Should differ by 1 at L
%figure,plot(numc)
reporttest('CELLCHUNK no overlap',allall(dnum==1))


[numc,latc,lonc]=cellchunk(num(3),lat(3),lon(3),200,'overlap');
numc(1,2:end)-numc(100,1:end-1);%Shoudl differ by 1 at L/2
%figure,plot(numc)
reporttest('CELLCHUNK with 50% overlap',allall(dnum==1))

