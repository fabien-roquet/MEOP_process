function [u_SA, u_CT, u_P] = gsw_internal_energy_first_derivatives_CT_exact(SA,CT,p)

% gsw_internal_energy_first_derivatives_CT_exact       first derivatives of
%                                       specific interal energy of seawater
%==========================================================================
%
% USAGE:
%  [u_SA, u_CT, u_P] = ...
%                   gsw_internal_energy_first_derivatives_CT_exact(SA,CT,p)
%
% DESCRIPTION:
%  Calculates the first order derivates of specific internal energy of 
%  seawater.
%
%  Note that this function uses the full Gibbs function.  There is an 
%  alternative to calling this function, namely 
%  gsw_internal_energy_first_derivatives(SA,CT,p), which uses the 
%  computationally efficient polynomial for specific volume in terms of 
%  SA, CT and p (Roquet et al., 2015).    
%
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  CT  =  Conservative Temperature                                [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & CT are MxN.
%
% OUTPUT:
%  u_SA = The first derivative of internal energy with respect to 
%           Absolute Salinity at constant CT & p.       
%                                          [ (J/kg)(g/kg)^-1 ] i.e. [ J/g ]
%  u_CT = The first derivative of internal energy with respect to 
%           Conservative Temperature at constant SA & p.    [ (J/kg) K^-1 ]
%  u_P = The first derivative of internal energy with respect to 
%           pressure at constant SA & CT.                  [ (J/kg) Pa^-1 ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker.                   [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05 (27th January 2015)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
    error('gsw_internal_energy_first_derivatives_CT_exact: requires three inputs')
end
if ~(nargout == 3)
    error('gsw_internal_energy_first_derivatives_CT_exact: requires three outputs')
end

[ms,ns] = size(SA);
[mt,nt] = size(CT); 
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_internal_energy_first_derivatives_CT_exact: SA and CT must have same dimensions')
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
    error('gsw_internal_energy_first_derivatives_CT_exact: Inputs array dimensions arguments do not agree')
end %if

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

db2Pa = 1e4;               %  dbar to Pa conversion factor 
P = (db2Pa.*p + gsw_P0);

[h_SA, h_CT] = gsw_enthalpy_first_derivatives_CT_exact(SA,CT,p);
v = gsw_specvol_CT_exact(SA,CT,p);
[v_SA, v_CT, v_P] = gsw_specvol_first_derivatives_CT_exact(SA,CT,p);

u_SA = h_SA - P.*v_SA;

u_CT = h_CT - P.*v_CT;

u_P = v - P.*v_P;

if transposed
    u_SA = u_SA.';
    u_CT = u_CT.';
    u_P = u_P.';
end

end
