function[ar]=latratio(lat,h)
%LATRATIO  Set plot aspect ratio for latitude / longitude plot.
%
%   LATRATIO(LAT) sets the aspect ratio of the current axis correctly a 
%   Cartesian plot centered about the latitude LAT, i.e. the x/y aspect 
%   ratio is set to COS(LAT).
%
%   Equal distances along the x- and y-axes then correspond to the same 
%   physical distance.
%
%   LAT is measured in degrees.
%
%   LATRATIO(LAT,H) does the same for the axis with handle H.
%
%   LATRATIO with no input arguments sets the aspect ratio of the current
%   axis using the midpoint of the y-axis limits.
%
%   AR=LATRATIO(LAT) alternately returns the correct aspect ratio in a 
%   two-element array AR, without modifying any plots. 
%
%   Usage: latratio(lat,h);
%          latratio(lat);
%          latratio;
%          ar=latration(lat);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2007--2013 J.M. Lilly --- type 'help jlab_license' for details
 
if nargin==0
   ax=axis;
   lat=frac(1,2)*(ax(3)+ax(4));
end

if nargin<2
    h=gca;
end

ar=[1./cos(jdeg2rad(lat)) 1 1];
if nargout==1
    ar=ar(1:2);
else
    set(h,'dataaspectratio',ar)
end
if nargout==0
    clear ar
end

