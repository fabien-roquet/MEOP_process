function[b]=maxmax(x)
%MAXMAX(X)=MAX(X(ISFINITE(X)))
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2012 J.M. Lilly --- type 'help jlab_license' for details        
b=max(x(isfinite(x(:))));
