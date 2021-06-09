function contours=list_contours(minZ,maxZ,Ncontours)
% function contours=list_contours(minZ,maxZ,Ncontours)
%
% determines a list of contours with nice roundings for contour plots.
% The number of contour selected should be around Ncontours.


dZ=(maxZ-minZ)/Ncontours;
n=0;
while dZ>=10, n=n+1; dZ=dZ/10; end
while dZ<1, n=n-1; dZ=dZ*10; end
dZ=floor(dZ)*10^n;
minZ=floor(minZ*10^-n)*10^n;
contours=minZ:dZ:maxZ;