function[c]=maternc(alpha)
%MATERNC  Returns the normalization function C_ALPHA for a Matern process.
%
%   MATERNC is a low-level function used in the JMATERN package.
%
%   MATERNC(ALPHA) returns (1/PI)*GAMMA(ALPHA-1/2)/GAMMA(ALPHA) where GAMMA
%   is the usual gamma function.
%
%   Usage: c=maternc(alpha);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details
 
c=frac(1,sqrt(pi))*frac(gamma(alpha-1/2),gamma(alpha));



