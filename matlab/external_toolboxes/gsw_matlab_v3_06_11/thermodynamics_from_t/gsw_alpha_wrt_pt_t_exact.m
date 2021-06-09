function alpha_wrt_pt_t_exact = gsw_alpha_wrt_pt_t_exact(SA,t,p)

% gsw_alpha_wrt_pt_t_exact                    thermal expansion coefficient
%                                     with respect to potential temperature
%==========================================================================
%
% USAGE:  
%  alpha_wrt_pt_t_exact = gsw_alpha_wrt_pt_t_exact(SA,t,p)
%
% DESCRIPTION:
%  Calculates the thermal expansion coefficient of seawater with respect to  
%  potential temperature, with a reference pressure of zero.
%   
% INPUT:
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA & t are MxN.
%
% OUTPUT:
%  alpha_wrt_pt_t_exact  =  thermal expansion coefficient           [ 1/K ]
%                           with respect to potential temperature,
%                           with a reference pressure of zero dbar.
%    
% AUTHOR: 
%  David Jackett, Trevor McDougall and Paul Barker     [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05 (27th January 2015)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See Eqn. (2.18.2) of this TEOS-10 manual.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3)
   error('gsw_alpha_wrt_pt_t_exact:  Requires three inputs')
end %if

[ms,ns] = size(SA);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt ~= ms | nt ~= ns)
    error('gsw_alpha_wrt_pt_t_exact: SA and t must have same dimensions')
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
    error('gsw_alpha_wrt_pt_t_exact: Inputs array dimensions arguments do not agree')
end %if

if ms == 1
    SA = SA.';
    t = t.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

pr0 = zeros(size(p)); 
pt0 = gsw_pt0_from_t(SA,t,p);

factor = gsw_gibbs(0,2,0,SA,pt0,pr0)./gsw_gibbs(0,2,0,SA,t,p);

alpha_wrt_pt_t_exact = factor.*(gsw_gibbs(0,1,1,SA,t,p)./gsw_gibbs(0,0,1,SA,t,p));

if transposed
    alpha_wrt_pt_t_exact = alpha_wrt_pt_t_exact.';
end

end
