function t_freezing = gsw_t_freezing_poly(SA,p,saturation_fraction)

% gsw_t_freezing_poly         in-situ temperature at which seawater freezes
%                                                                    (poly)
%==========================================================================
%
% USAGE:
%  t_freezing = gsw_t_freezing_poly(SA,p,saturation_fraction)
%
% DESCRIPTION:
%  Calculates the in-situ temperature at which seawater freezes from a 
%  comptationally efficient polynomial.
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
% OPTIONAL:
%  saturation_fraction = the saturation fraction of dissolved air in 
%                        seawater
%  (i.e., saturation_fraction must be between 0 and 1, and the default 
%    is 0, air free) 
%
%  p & saturation_fraction (if provided) may have dimensions 1x1 or Mx1 or 
%  1xN or MxN, where SA is MxN.
%
% OUTPUT:
%  t_freezing = in-situ temperature at which seawater freezes.    [ deg C ]
%               (ITS-90)                
%
% AUTHOR: 
%  Trevor McDougall, Paul Barker and Rainer Feistal    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05 (27th January 2015)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org.
%    See sections 3.33 and 3.34 of this TEOS-10 Manual.  
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014: 
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation. 
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2 | nargin == 3) 
   error('gsw_t_freezing_poly: Requires either two or three inputs')
end %if

if ~exist('saturation_fraction','var')
    saturation_fraction = 0;
end
    
if (saturation_fraction < 0 | saturation_fraction > 1)
   error('gsw_t_freezing_poly: saturation_fraction MUST be between zero and one.')
end
    
[ms,ns] = size(SA);
[mp,np] = size(p);
[msf,nsf] = size(saturation_fraction);

if (mp == 1) & (np == 1)                    % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)                            % p is row vector,
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (np == 1)                         % p is column vector,
    p = p(:,ones(1,ns));                            % copy across each row.
elseif (ns == mp) & (np == 1)               % p is a transposed row vector,
    p = p.';                                               % transposed then
    p = p(ones(1,ms), :);                          % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_t_freezing_poly: Inputs array dimensions arguments do not agree')
end

if (msf == 1) & (nsf == 1)                                    % saturation_fraction scalar
    saturation_fraction = saturation_fraction*ones(size(SA));         % fill to size of SA
elseif (ns == nsf) & (msf == 1)                        % saturation_fraction is row vector,
    saturation_fraction = saturation_fraction(ones(1,ms), :);      % copy down each column.
elseif (ms == msf) & (nsf == 1)                     % saturation_fraction is column vector,
    saturation_fraction = saturation_fraction(:,ones(1,ns));        % copy across each row.
elseif (ns == msf) & (nsf == 1)           % saturation_fraction is a transposed row vector,
    saturation_fraction = saturation_fraction.';                           % transposed then
    saturation_fraction = saturation_fraction(ones(1,ms), :);      % copy down each column.
elseif (ms == msf) & (ns == nsf)
    % ok
else
    error('gsw_t_freezing_poly: Inputs array dimensions arguments do not agree')
end

if ms == 1
    SA = SA.';
    p = p.';
    saturation_fraction = saturation_fraction.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SA(SA < 0) = 0; % This line ensure that SA is non-negative.

CT_freezing = gsw_CT_freezing_poly(SA,p,saturation_fraction);
t_freezing = gsw_t_from_CT(SA,CT_freezing,p);

%--------------------------------------------------------------------------
% This function, gsw_t_freezing_poly, calculates the in-situ freezing 
% temperature, t_freezing, of seawater by first evaluating a polynomial of 
% the Conservative Temperature at which seawater freezes, CT_freezing, 
% using the GSW function gsw_CT_freezing_poly.  The in-situ freezing temperature
% is then calculated using the GSW function gsw_t_from_CT.  However, if one
% wanted to compute the in-situ freezing temperature directly from a single
% polynomial expression without first calculating the Conservative 
% Temperature at the freezing point, the following lines of code achieve 
% this.  The error of the following fit is similar to that of the present 
% function, gsw_t_freezing_poly, and ranges between -8e-4 K and 3e-4 K when 
% compared with the in-situ freezing temperature evaluated by Newton-
% Raphson iteration of the equality of the chemical potentials of water in 
% seawater and in ice.  (Note that the in-situ freezing temperature can be 
% found by this exact method using the function gsw_t_freezing).  
% 
% c0 = 0.002519;
% 
% c1 = -5.946302841607319;
% c2 =  4.136051661346983;
% c3 = -1.115150523403847e1;
% c4 =  1.476878746184548e1;
% c5 = -1.088873263630961e1;
% c6 =  2.961018839640730;
%     
% c7 = -7.433320943962606;
% c8 = -1.561578562479883;
% c9 =  4.073774363480365e-2;
% 
% c10 =  1.158414435887717e-2;
% c11 = -4.122639292422863e-1;
% c12 = -1.123186915628260e-1;
% c13 =  5.715012685553502e-1;
% c14 =  2.021682115652684e-1;
% c15 =  4.140574258089767e-2;
% c16 = -6.034228641903586e-1;
% c17 = -1.205825928146808e-2;
% c18 = -2.812172968619369e-1;
% c19 =  1.877244474023750e-2;
% c20 = -1.204395563789007e-1;
% c21 =  2.349147739749606e-1;
% c22 =  2.748444541144219e-3;
% 
% SA_r = SA.*1e-2;
% x = sqrt(SA_r);
% p_r = p.*1e-4;
% 
% t_freezing = c0 ...
%     + SA_r.*(c1 + x.*(c2 + x.*(c3 + x.*(c4 + x.*(c5 + c6.*x))))) ...
%     + p_r.*(c7 + p_r.*(c8 + c9.*p_r)) ...
%     + SA_r.*p_r.*(c10 + p_r.*(c12 + p_r.*(c15 + c21.*SA_r)) + SA_r.*(c13 + c17.*p_r + c19.*SA_r) ...
%     + x.*(c11 + p_r.*(c14 + c18.*p_r)  + SA_r.*(c16 + c20.*p_r + c22.*SA_r)));
% 
% Adjust for the effects of dissolved air
% t_freezing = t_freezing  - saturation_fraction.*(1e-3).*(2.4 - SA./70.33008);
% 
%---------------This is the end of the alternative code--------------------

t_freezing(p > 10000 | SA > 120 | ...
    p + SA.*71.428571428571402 > 13571.42857142857) = NaN;

if transposed
    t_freezing = t_freezing.';
end

end