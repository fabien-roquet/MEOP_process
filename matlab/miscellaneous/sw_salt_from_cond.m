
function S = sw_salt_from_cond(C,T,P)

% SW_SALT    Salinity from C, T, P
%=========================================================================
% SW_SALT  $Id: sw_salt.m,v 1.2 2004/03/18 20:27:03 pen078 Exp $
%          Copyright (C) CSIRO, Phil Morgan 1993.
%
% USAGE: S = sw_salt_from_cond(C,T,P)
%
% DESCRIPTION:
%   Calculates Salinity from Conductivity. UNESCO 1983 polynomial.
%
% INPUT:
%   C    = Conductivity [mS/cm]
%   T    = temperature  [degree C (ITS-90)]
%   P    = pressure     [db]
%
% OUTPUT:
%   S    = salinity    [psu      (PSS-78)]
%
% AUTHOR:  Phil Morgan 93-04-17, Lindsay Pender (Lindsay.Pender@csiro.au)
%
% DISCLAIMER:
%   This software is provided "as is" without warranty of any kind.
%   See the file sw_copy.m for conditions of use and licence.
%
% REFERENCES:
%    Fofonoff, P. and Millard, R.C. Jr
%    Unesco 1983. Algorithms for computation of fundamental properties of
%    seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
%=========================================================================

% Modifications
% 03-12-12. Lindsay Pender, Converted to ITS-90.
% 20-04-2015. Fabien Roquet.

% CALLER: general purpose
% CALLEE: sw_sals.m sw_salrt.m sw_salrp.m


%----------------------------------
% CHECK INPUTS ARE SAME DIMENSIONS
%----------------------------------
[mc,nc] = size(C);
[mt,nt] = size(T);
[mp,np] = size(P);

if ~(mc==mt | mc==mp | nc==nt | nc==np)
  error('sw_salt_from_cond.m: C,T,P must all have the same dimensions')
end %if

%-------
% BEGIN
%-------
R  = C./sw_c3515;
rt = sw_salrt(T);
Rp = sw_salrp(R,T,P);
Rt = R./(Rp.*rt);
S  = sw_sals(Rt,T);

