function specvol_ice = gsw_specvol_from_pot_enthalpy_ice(pot_enthalpy_ice,p)

% gsw_specvol_from_pot_enthalpy_ice                    specific volume from
%                                                 potential enthalpy of ice                              
%==========================================================================
%
% USAGE:
%  specvol_ice = gsw_specvol_from_pot_enthalpy_ice(pot_enthalpy_ice,p)
%
% DESCRIPTION:
%  Calculates the specific volume of ice from the potential enthalpy 
%  of ice.  The reference sea pressure of the potential enthalpy is zero
%  dbar. 
%
% INPUT:
%  pot_enthalpy_ice  =  potential enthalpy of ice                  [ J/kg ]
%  p  =  sea pressure                                              [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar ) 
%
% OUTPUT:
%  specvol_ice  =  specific volume                               [ m^3/kg ]
%
% AUTHOR: 
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05 (27th January 2015)
%
% REFERENCES:
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%
%  McDougall, T.J., P.M. Barker, R. Feistel and B.K. Galton-Fenzi, 2014: 
%   Melting of Ice and Sea Ice into Seawater and Frazil Ice Formation. 
%   Journal of Physical Oceanography, 44, 1751-1775.
%
%  McDougall T. J. and S. J. Wotherspoon, 2014: A simple modification of 
%   Newton's method to achieve convergence of order 1 + sqrt(2).  Applied 
%   Mathematics Letters, 29, 20-25.  
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 2)
   error('gsw_specvol_from_pot_enthalpy_ice:  Requires two inputs')
end

[mh,nh] = size(pot_enthalpy_ice);
[mp,np] = size(p);

if (mp == 1) & (np == 1)           
    p = p*ones(mh,nh);
elseif (nh == np) & (mp == 1)         
    p = p(ones(1,mh), :);             
elseif (mh == mp) & (np == 1)        
    p = p(:,ones(1,nh));               
elseif (nh == mp) & (np == 1)        
    p = p.';                     
    p = p(ones(1,mh), :);          
elseif (mh == mp) & (nh == np)
    % ok
else
    error('gsw_specvol_from_pot_enthalpy_ice: Inputs array dimensions arguments do not agree')
end

if mh == 1
    pot_enthalpy_ice = pot_enthalpy_ice.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

pt0_ice = gsw_pt_from_pot_enthalpy_ice(pot_enthalpy_ice);

t_ice = gsw_t_from_pt0_ice(pt0_ice,p);

specvol_ice = gsw_specvol_ice(t_ice,p);

if transposed
   specvol_ice = specvol_ice.';
end

end
