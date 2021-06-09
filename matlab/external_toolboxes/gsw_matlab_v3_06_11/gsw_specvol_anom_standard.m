function specvol_anom = gsw_specvol_anom_standard(SA,CT,p)

% gsw_specvol_anom_standard                    specific volume anomaly with
%                         reference of SA = SSO & CT = 0 (75-term equation)
%==========================================================================
%
% USAGE:
%  specvol_anom = gsw_specvol_anom_standard(SA,CT,p)
%
% DESCRIPTION:
%  Calculates specific volume anomaly from Absolute Salinity, Conservative
%  Temperature and pressure. It uses the computationally-efficient
%  expression for specific volume as a function of SA, CT and p (Roquet
%  et al., 2015).  The reference value to which the anomally is calculated
%  has an Absolute Salinity of SSO and Conservative Temperature equal to
%  0 degress C.
%
%  Note that this 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:
%  SA     =  Absolute Salinity                                     [ g/kg ]
%  CT     =  Conservative Temperature (ITS-90)                    [ deg C ]
%  p      =  sea pressure                                          [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & CT  need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  specvol_anom  =  specific volume anomaly                      [ m^3/kg ]
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05 (27th January 2015)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of
%   seawater - 2010: Calculation and use of thermodynamic properties.
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (3.7.3) of this TEOS-10 Manual.
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling.
%
% The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
    error('gsw_specvol_anom_standard: Requires three inputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_specvol_anom_standard: SA and CT must have same dimensions')
end

if (mp == 1) & (np == 1)              % p scalar - fill to size of SA
    p = p*ones(size(SA));
elseif (ns == np) & (mp == 1)         % p is row vector,
    p = p(ones(1,ms), :);              % copy down each column.
elseif (ms == mp) & (np == 1)         % p is column vector,
    p = p(:,ones(1,ns));               % copy across each row.
elseif (ns == mp) & (np == 1)          % p is a transposed row vector,
    p = p.';                              % transposed then
    p = p(ones(1,ms), :);                % copy down each column.
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_specvol_anom_standard: Inputs array dimensions arguments do not agree')
end

if ms == 1
    SA = SA.';
    CT = CT.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

% This line ensures that SA is non-negative.
SA(SA < 0) = 0;

%deltaS = 24;
sfac = 0.0248826675584615;                   % sfac = 1/(40*(35.16504/35)).
offset = 5.971840214030754e-1;                      % offset = deltaS*sfac.

x2 = sfac.*SA;
xs = sqrt(x2 + offset);
ys = CT.*0.025;
z = p.*1e-4;

%     v000 =  1.0769995862e-3;
%     v001 = -6.0799143809e-5;
%     v002 =  9.9856169219e-6;
%     v003 = -1.1309361437e-6;
%     v004 =  1.0531153080e-7;
%     v005 = -1.2647261286e-8;
%     v006 =  1.9613503930e-9;
v010 = -1.5649734675e-5;
v011 =  1.8505765429e-5;
v012 = -1.1736386731e-6;
v013 = -3.6527006553e-7;
v014 =  3.1454099902e-7;
v020 =  2.7762106484e-5;
v021 = -1.1716606853e-5;
v022 =  2.1305028740e-6;
v023 =  2.8695905159e-7;
v030 = -1.6521159259e-5;
v031 =  7.9279656173e-6;
v032 = -4.6132540037e-7;
v040 =  6.9111322702e-6;
v041 = -3.4102187482e-6;
v042 = -6.3352916514e-8;
v050 = -8.0539615540e-7;
v051 =  5.0736766814e-7;
v060 =  2.0543094268e-7;
v100 = -3.1038981976e-4;
v101 =  2.4262468747e-5;
v102 = -5.8484432984e-7;
v103 =  3.6310188515e-7;
v104 = -1.1147125423e-7;
v110 =  3.5009599764e-5;
v111 = -9.5677088156e-6;
v112 = -5.5699154557e-6;
v113 = -2.7295696237e-7;
v120 = -3.7435842344e-5;
v121 = -2.3678308361e-7;
v122 =  3.9137387080e-7;
v130 =  2.4141479483e-5;
v131 = -3.4558773655e-6;
v132 =  7.7618888092e-9;
v140 = -8.7595873154e-6;
v141 =  1.2956717783e-6;
v150 = -3.3052758900e-7;
v200 =  6.6928067038e-4;
v201 = -3.4792460974e-5;
v202 = -4.8122251597e-6;
v203 =  1.6746303780e-8;
v210 = -4.3592678561e-5;
v211 =  1.1100834765e-5;
v212 =  5.4620748834e-6;
v220 =  3.5907822760e-5;
v221 =  2.9283346295e-6;
v222 = -6.5731104067e-7;
v230 = -1.4353633048e-5;
v231 =  3.1655306078e-7;
v240 =  4.3703680598e-6;
v300 = -8.5047933937e-4;
v301 =  3.7470777305e-5;
v302 =  4.9263106998e-6;
v310 =  3.4532461828e-5;
v311 = -9.8447117844e-6;
v312 = -1.3544185627e-6;
v320 = -1.8698584187e-5;
v321 = -4.8826139200e-7;
v330 =  2.2863324556e-6;
v400 =  5.8086069943e-4;
v401 = -1.7322218612e-5;
v402 = -1.7811974727e-6;
v410 = -1.1959409788e-5;
v411 =  2.5909225260e-6;
v420 =  3.8595339244e-6;
v500 = -2.1092370507e-4;
v501 =  3.0927427253e-6;
v510 =  1.3864594581e-6;
v600 =  3.1932457305e-5;

xy_part_0_SSO_0 = -1.043382007156129e-4;  % xy_part_0 evaluated at SA=SSO and CT=0.
xy_part_1_SSO_0 =  1.574001169739070e-5;  % xy_part_1 evaluated at SA=SSO and CT=0.
xy_part_2_SSO_0 = -2.854887955972872e-6;  % xy_part_2 evaluated at SA=SSO and CT=0.
xy_part_3_SSO_0 =  4.652181957231689e-7;  % xy_part_3 evaluated at SA=SSO and CT=0.
xy_part_4_SSO_0 = -1.352520752723288e-7;  % xy_part_4 evaluated at SA=SSO and CT=0.

xy_part_0 = xs.*(v100 + xs.*(v200 + xs.*(v300 + xs.*(v400 + xs.*(v500 ...
    + v600.*xs))))) + ys.*(v010 + xs.*(v110 + xs.*(v210 + xs.*(v310 + xs.*(v410 ...
    + v510.*xs)))) + ys.*(v020 + xs.*(v120 + xs.*(v220 + xs.*(v320 + v420.*xs))) ...
    + ys.*(v030 + xs.*(v130 + xs.*(v230 + v330.*xs)) + ys.*(v040 + xs.*(v140 ...
    + v240*xs) + ys.*(v050 + v150.*xs + v060.*ys)))));

xy_part_1 = xs.*(v101 + xs.*(v201 + xs.*(v301 + xs.*(v401 + v501.*xs)))) ...
    + ys.*(v011 + xs.*(v111 + xs.*(v211 + xs.*(v311 + v411.*xs))) + ys.*(v021 ...
    + xs.*(v121 + xs.*(v221 + v321.*xs)) + ys.*(v031 + xs.*(v131 + v231.*xs) ...
    + ys.*(v041 + v141.*xs + v051.*ys))));

xy_part_2 = xs.*(v102 + xs.*(v202 + xs.*(v302 + v402.*xs))) ...
    + ys.*(v012 + xs.*(v112 + xs.*(v212 + v312.*xs)) + ys.*(v022 + xs.*(v122 ...
    + v222.*xs) + ys.*(v032 + v132.*xs + v042.*ys)));

xy_part_3 = xs.*(v103 + v203.*xs) + ys.*(v013 + v113.*xs + v023.*ys);

xy_part_4 = v104.*xs + v014.*ys;

specvol_anom = xy_part_0 - xy_part_0_SSO_0 + z.*(xy_part_1 - xy_part_1_SSO_0 ...
    + z.*(xy_part_2 - xy_part_2_SSO_0 + z.*(xy_part_3 - xy_part_3_SSO_0 ...
    + z.*(xy_part_4 - xy_part_4_SSO_0))));

% Note that in the above xy_part_(0 to 4) the constant terms (v000 to v004)
% are not included and the pressure terms that are to the 5th and 6th power
% are not evaluated as they cancel.

%--------------------------------------------------------------------------
% This function calculates specvol_anom using the computationally
% efficient expression for specific volume in terms of SA, CT and p.  If
% one wanted to compute specvol_anom from SA, CT, and p with the full
% TEOS-10 Gibbs function, the following lines of code will enable this.
%
%    t = gsw_t_from_CT(SA,CT,p);
%    specvol_anom = gsw_specvol_anom_standard_t_exact(SA,t,p);
%
%    or call the following, it is identical to the lines above.
%
%    specvol_anom = gsw_specvol_anom_standard_CT_exact(SA,CT,p)
%
%-----------------This is the end of the alternative code------------------

if transposed
    specvol_anom = specvol_anom.';
end

end
