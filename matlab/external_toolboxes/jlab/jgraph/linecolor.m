function[varargout]=linecolor(varargin)
%LINECOLOR  Set line colors based on a property value within a colormap.
%
%   LINECOLOR(H,C) set the lines with handles H to the colors C, determined
%   by looking up the values of C within the 'jet' colormap.  
%
%   H and C should be arrays of the same size.
%
%   LINECOLOR(H,C,CMIN,CMAX) uses the values CMIN and CMAX as the lower and
%   uppermost values in the colormap, respectively.  The default behavior 
%   is to set CMIN and CMAX to the minimum and maximum values within C.
%
%   LINECOLOR(...,MAP) alternately uses the colormap with the name MAP.
%
%   'linecolor --t' runs a test.
%
%   Usage: linecolor(h,c);
%          linecolor(h,c,map);
%          linecolor(h,c,cmin,cmax);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013 J.M. Lilly --- type 'help jlab_license' for details

h=varargin{1};
c=varargin{2};
if ischar(varargin{end})
    str=varargin{end};
    varargin=varargin(1:end-1);
else
    str='jet';
end
varargin=varargin(3:end);

if length(varargin)==2
    cmin=varargin{1};
    cmax=varargin{2};
else
    cmin=minmin(c);
    cmax=maxmax(c);
end

m=frac(64-1,cmax-cmin);


b=1-m*cmin;
y=ceil(m*c+b);

map=colormap(str);

y(y>64)=64;
y(y<1)=1;

set(h,'visible','off')
for i=1:length(h),
    set(h(i),'color',map(y(i),:));
end

set(h,'visible','on')

