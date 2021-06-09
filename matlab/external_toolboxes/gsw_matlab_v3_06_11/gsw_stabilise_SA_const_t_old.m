function SA_out = gsw_stabilise_SA_const_t(SA_in,t,p,opt_1,opt_2)

% gsw_stabilise_SA_const_t             adjusts salinities (SA) to produce a
%                             water column stablised to be neutral, in-situ
%                           temperature remains constant (75-term equation)
%==========================================================================
%
% USAGE:
%  SA_out = gsw_stabilise_SA_const_t(SA_in,t,p,{opt_1,opt_2})
%
% DESCRIPTION:
%  This function stabilises a water column. This is achieved by minimally 
%  adjusting only the Absolute Salinity SA values such that the minimum 
%  stability is made to be within at least 1 x 10^-9 s^-2 of the desired 
%  minimum Nsquared min_Nsquared, the default value is which is about 1/5th 
%  of the square of earth's rotation rate. There are no changes made to
%  either in-situ temperature or pressure.
%
%  This programme requires either the Optimization toolbox or Tomlab CPLEX.
%  if there are a up to several hundred data points in the cast then 
%  Matlab's Optimization toolbox produces reasonable results, but if there 
%  are thousands of bottles in the cast or the best possible output is  
%  wanted then the CPLEX solver is required. This programme will determine
%  if Tomlab or the Optimization toolbox is available to the user, if both
%  are available it will use Tomlab.
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA_in  =  uncorrected Absolute Salinity                         [ g/kg ]
%  t      =  in-situ temperature (ITS-90)                         [ deg C ]
%  p      =  sea pressure                                          [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
% OPTIONAL:
%  opt_1 = Nsquared_lowerlimit                                    [ 1/s^2 ]
%  Note. If Nsquared_lowerlimit is not supplied, a default minimum 
%   stability of 1 x 10^-9 s^-2 will be applied.
%  or,
%  opt_1 =  longitude in decimal degrees                     [ 0 ... +360 ]
%                                                     or  [ -180 ... +180 ]
%  opt_2 =  latitude in decimal degrees north               [ -90 ... +90 ] 
%
%  SA_in & t need to have the same dimensions.
%  p may have dimensions 1x1 or Mx1 or 1xN or MxN, where SA_in & t are MxN.
%  opt_1 equal to Nsquared_lowerlimit, if provided, may have dimensions 1x1 
%  or (M-1)x1 or 1xN or (M-1)xN, where SA_in & t are MxN.
%  opt_1 equal to long & opt_2 equal to lat, if provided, may have
%  dimensions 1x1 or (M-1)x1 or 1xN or (M-1)xN, where SA_in & t are MxN.
%
% OUTPUT:
%  SA_out =  corrected stabilised Absolute Salinity                [ g/kg ]
%
% AUTHOR:
%  Paul Barker and Trevor McDougall                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05.5 (3rd June, 2016)
%
% REFERENCES:
%  Barker, P.M., and T.J. McDougall, 2016: Stabilisation of hydrographic 
%    profiles.  J. Atmosph. Ocean. Tech., submitted.
%
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org 
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling, 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%  The Tomlab software is available from http://www.tomopt.com
%
%==========================================================================

%--------------------------------------------------------------------------
% Check if necessary software exists
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 3 || nargin == 4 || nargin == 5)
    error('gsw_stabilise_SA_const_t:  Requires three or four or five inputs')
end
     
[ms,ns] = size(SA_in);
[mt,nt] = size(t);
[mp,np] = size(p);

if (mt ~= ms || nt ~= ns)
    error('gsw_stabilise_SA_const_t: SA_in and t must have same dimensions')
end
if (ms*ns == 1)
    error('gsw_stabilise_SA_const_t:  There must be at least 3 bottles')
end

if (mp == 1) & (np == 1)
    error('gsw_stabilise_SA_const_t:  There must be at least 3 bottles')
elseif (ns == np) & (mp == 1)
    p = p(ones(1,ms), :);
elseif (ms == mp) & (np == 1)
    p = p(:,ones(1,ns));
elseif (ns == mp) & (np == 1)
    p = p.'; 
    p = p(ones(1,ms), :);
elseif (ms == np) & (mp == 1)
     p = p.';  
     p = p(:,ones(1,ns));
elseif (ms == np) & (ns == mp)
     p = p.';   
elseif (ms == mp) & (ns == np)
    % ok
else
    error('gsw_stabilise_SA_const_t: Inputs array dimensions arguments do not agree')
end 

if ms == 1
    SA_in = SA_in.';
    t = t.';
    p = p.';
    transposed = 1;
else
    transposed = 0;
end

[mp,number_profiles] = size(p); 

if nargin == 4
    Nsquared_lowerlimit = opt_1;
    if transposed
        Nsquared_lowerlimit = Nsquared_lowerlimit.';
    end
    [mN2,nN2] = size(Nsquared_lowerlimit);
    if (mN2 == 1) & (nN2 == 1)              
        Nsquared_lowerlimit_tmp = Nsquared_lowerlimit*ones(mp,number_profiles);
    elseif (number_profiles == nN2) & (mN2 == 1)     
        Nsquared_lowerlimit_tmp = Nsquared_lowerlimit(ones(1,mp),:);        
    elseif (mp == (mN2+1)) & (nN2 == 1)        
        Nsquared_lowerlimit_tmp = NaN(mp,number_profiles);
        Nsquared_lowerlimit_tmp(2:end,:) = Nsquared_lowerlimit(:,ones(1,number_profiles)); 
    elseif (mp == (mN2+1)) & (number_profiles == nN2)
        Nsquared_lowerlimit_tmp = NaN(mp,number_profiles);
        Nsquared_lowerlimit_tmp(2:end,:) = Nsquared_lowerlimit;
    else
        error('gsw_stabilise_SA_const_t: Inputs array dimensions arguments do not agree')
    end
end

if nargin == 5
    long = opt_1;
    lat = opt_2;
    
    [mlo,nlo] = size(long);
    long(long < 0) = long(long < 0) + 360;
    
    if (mlo == 1) & (nlo == 1) 
        long = long*ones(mp,number_profiles);
    elseif (number_profiles == nlo) & (mlo == 1)
        long = long(ones(1,mp), :); 
    elseif (mp == mlo) & (nlo == 1)
        long = long(:,ones(1,number_profiles));
    elseif (number_profiles == mlo) & (nlo == 1)
        long = long.'; 
        long = long(ones(1,mp), :); 
    elseif (mp == nlo) & (mlo == 1)
        long = long.'; 
        long = long(:,ones(1,number_profiles));
    elseif (mp == mlo) & (number_profiles == nlo)
        % ok
    else
        error('gsw_stabilise_SA_const_t: Inputs array dimensions arguments do not agree')
    end
    
    [mla,nla] = size(lat);
    
    if (mla == 1) & (nla == 1)
        lat = lat*ones(mp,number_profiles);
    elseif (number_profiles == nla) & (mla == 1) 
        lat = lat(ones(1,mp), :);
    elseif (mp == mla) & (nla == 1) 
        lat = lat(:,ones(1,number_profiles));
    elseif (number_profiles == mla) & (nla == 1) 
        lat = lat.'; 
        lat = lat(ones(1,mp), :);
    elseif (mp == mla) & (number_profiles == nla)
        % ok
    else
        error('gsw_stabilise_SA_const_t: Inputs array dimensions arguments do not agree')
    end
    
    Nsquared_lowerlimit_tmp = gsw_Nsquared_lowerlimit(p,long,lat);

    clear long mla nla lat mlo mlo 
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

Nsquared_lowerlimit_default = 1e-9;

% db2Pa = 1e4;
% grav = 9.7963 (Griffies, 2004)
c =1.250423402612047e+02; % c = 1.2*db2Pa./(grav.^2);
rho0 = 1e3;
A_0 = eye(mp) - diag(ones(mp-1,1),1); A_0(mp,:) = [];

%--------------------------------------------------------------------------
% set TEOS-10 limits
%--------------------------------------------------------------------------
t(p < 100 & (t > 80 | t < -12)) = NaN;
t(p >= 100 & (t > 40 | t < -12)) = NaN;
t(SA_in > 120 | t < -12 | t > 80 | p > 12000) = NaN;
t(p < -1.5 | p > 12000) = NaN;
%--------------------------------------------------------------------------

SA_out = NaN(mp,number_profiles);
opts = optimoptions('quadprog','Algorithm','interior-point-convex','Display','off');

for Iprofile = 1:number_profiles
       
    [Inn] = find(~isnan(SA_in(:,Iprofile) + t(:,Iprofile) + p(:,Iprofile)));
    
    if length(Inn) < 2
        SA_out(Inn,Iprofile) = SA_in(Inn,Iprofile);
        continue
    end
    
    SA_tmp = SA_in(Inn,Iprofile);
    t_tmp = t(Inn,Iprofile);
    p_tmp = p(Inn,Iprofile);
    
    Ishallow = 1:(length(Inn)-1);
    Ideep = 2:length(Inn);
    d_p = (p_tmp(Ideep) - p_tmp(Ishallow));
    
    if any(d_p <= 0)
        warning('gsw_stabilise_SA_const_t: pressure must be monotonic')
        continue
    end
    
    % calculate the minimum Nsquared
    [N2_tmp,N2_alpha_tmp,N2_beta_tmp,dSA_tmp,dCT_tmp,dp_tmp] ...
        = gsw_Nsquared_min_const_t(SA_tmp,t_tmp,p_tmp);
    
    pl = length(p_tmp);
    
    %--------------------------------------------------------------------------
    % Set the Nsquared lower limit
    %--------------------------------------------------------------------------
    if ~exist('Nsquared_lowerlimit_tmp','var')
        Nsquared_lowerlimit = Nsquared_lowerlimit_default*ones(pl-1,1); %default
    else
        dummy = squeeze(Nsquared_lowerlimit_tmp(Inn,Iprofile));
        Nsquared_lowerlimit = dummy(2:end);
    end
    %--------------------------------------------------------------------------
    
    [Iunstable] = find((N2_tmp - Nsquared_lowerlimit) < 0);
    
    if isempty(Iunstable)
        SA_out(:,Iprofile) = SA_in(:,Iprofile);
        continue;
    end
    
    Name = 'stabilise the water column by adjusting SA while keeping t constant';
    
    b_U = dSA_tmp - (N2_alpha_tmp./N2_beta_tmp).*dCT_tmp ...
        - c/rho0*(Nsquared_lowerlimit.*dp_tmp./N2_beta_tmp) ;
    % Note that c = 1.2*db2Pa./(grav.^2);
    
    % Matlab solver
    H = speye(pl);
    A = A_0(1:pl-1,1:pl);
    x_L = - SA_tmp;
    x = quadprog(H, [], A, b_U, [], [], x_L, [], [], opts);
    
    SA_tmp = SA_tmp + x;
    
    %--------------------------------------------------------------------------
    
    
    SA_out(Inn,Iprofile) = SA_tmp;
    
end

if transposed
    SA_out = SA_out.';
end

end


function [N2, N2_alpha, N2_beta, dSA, dCT, dp] = gsw_Nsquared_min_const_t(SA,t,p)

% gsw_Nsquared_min_const_t               buoyancy (Brunt-Vaisala) frequency 
%                 squared (N^2) from in-situ temperature (75-term equation)
%==========================================================================
% 
% USAGE:  
%  [N2, N2_p, N2_specvol, N2_alpha, N2_beta, dSA, dCT, dp, N2_beta_ratio] = ...
%                                    gsw_Nsquared_min_const_t(SA,t,p,{lat})
%
% DESCRIPTION:
%  Calculates the minimum buoyancy frequency squared (N^2)(i.e. the 
%  Brunt-Vaisala frequency squared) between two bottles from the equation,
%
%           2      2    beta.dSA - alpha.dCT
%         N   =  g  . -------------------------
%                         specvol_local.dP
%
%  The pressure increment, dP, in the above formula is in Pa, so that it is
%  10^4 times the pressure increment dp in dbar. 
%
%  Note that this 75-term equation has been fitted in a restricted range of 
%  parameter space, and is most accurate inside the "oceanographic funnel" 
%  described in McDougall et al. (2003).  The GSW library function 
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if 
%  some of one's data lies outside this "funnel".  
%
% INPUT:  
%  SA  =  Absolute Salinity                                        [ g/kg ]
%  t   =  in-situ temperature (ITS-90)                            [ deg C ]
%  p   =  sea pressure                                             [ dbar ]
%         ( i.e. absolute pressure - 10.1325 dbar )
%
%  SA, t, p need to have the same dimensions. 
%
% OUTPUT:
%  N2         =  minimum Brunt-Vaisala Frequency squared          [ 1/s^2 ]
%  N2_alpha   =  thermal expansion coefficient with respect         [ 1/K ]
%                to Conservative Temperature at the minimum N2
%  N2_beta    =  saline contraction coefficient at constant        [ kg/g ]
%                Conservative Temperature at the minimum N2
%  dSA        =  salinity difference between bottles               [ g/kg ]
%  dCT        =  Conservative Temperature difference between      [ deg C ]
%                bottles
%  dp         =  pressure difference between bottles               [ dbar ]
%
% AUTHOR:  
%  Trevor McDougall and Paul Barker                    [ help@teos-10.org ]
%
% VERSION NUMBER: 3.05.5 (6th June, 2016)
%
% REFERENCES:
%  Griffies, S. M., 2004: Fundamentals of Ocean Climate Models. Princeton, 
%   NJ: Princeton University Press, 518 pp + xxxiv.
%   
%  IOC, SCOR and IAPSO, 2010: The international thermodynamic equation of 
%   seawater - 2010: Calculation and use of thermodynamic properties.  
%   Intergovernmental Oceanographic Commission, Manuals and Guides No. 56,
%   UNESCO (English), 196 pp.  Available from http://www.TEOS-10.org
%    See section 3.10 and Eqn. (3.10.2) of this TEOS-10 Manual. 
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003: 
%   Accurate and computationally efficient algorithms for potential 
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling, 90, pp. 29-43.
%
%   The software is available from http://www.TEOS-10.org
%
%==========================================================================

[mp,number_profiles] = size(p); 
grav2 = 9.7963^2;

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

db2Pa = 1e4;
rho0 = 1e3;
CT = t;

Ishallow = 1:(mp-1);
Ideep = 2:mp;

dSA = SA(Ideep,:) - SA(Ishallow,:);
dCT = CT(Ideep,:) - CT(Ishallow,:);
dp = p(Ideep,:) - p(Ishallow,:);

[specvol_bottle alpha_bottle beta_bottle]= gsw_specvol_alpha_beta(SA,CT,p);

N2_shallow = (rho0*grav2./(db2Pa.*dp)).*(beta_bottle(Ishallow,:).*dSA - alpha_bottle(Ishallow,:).*dCT);
N2_deep = (rho0*grav2./(db2Pa.*dp)).*(beta_bottle(Ideep,:).*dSA - alpha_bottle(Ideep,:).*dCT);

N2 = nan(mp-1,number_profiles);
N2_alpha = nan(mp-1,number_profiles);
N2_beta = nan(mp-1,number_profiles);

for Iprofile = 1:number_profiles
    dummy_N2 = [N2_shallow(:,Iprofile),N2_deep(:,Iprofile)];
    dummy_alpha = [alpha_bottle(Ishallow,Iprofile),alpha_bottle(Ideep,Iprofile)];
    dummy_beta = [beta_bottle(Ishallow,Iprofile),beta_bottle(Ideep,Iprofile)];

    [N2(:,Iprofile),IN2] = min(dummy_N2,[],2);
    
    for Ibottle = 1:mp-1
        N2_alpha(Ibottle,Iprofile) = dummy_alpha(Ibottle,IN2(Ibottle));
        N2_beta(Ibottle,Iprofile) = dummy_beta(Ibottle,IN2(Ibottle));
    end
end

end
