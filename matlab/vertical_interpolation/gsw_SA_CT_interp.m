function [SA_i, CT_i] = gsw_SA_CT_interp(SA,CT,p,p_i,factor)

% gsw_SA_CT_interp                                  SA and CT interpolation
%                                                          to p_i on a cast
%==========================================================================
%
% USAGE:
%  [SA_i, CT_i] = gsw_SA_CT_interp(SA,CT,p,p_i)
%
% DESCRIPTION:
%  Interpolate Absolute Salinity and Conservative Temperature values to
%  arbitrary pressures using the SA-CT diagram.  Any interpolated bottles 
%  that have pressures shallower than the shallowest observed bottle are
%  set equal to the shallowest observed bottle.
%
%  Note that this interpolation scheme requires at least four observed
%  bottles on the cast.
%
%  This programme is based on the computationally-efficient expression
%  for specific volume in terms of SA, CT and p (Roquet et al., 2015).
%
%  Note that this 75-term equation has been fitted in a restricted range of
%  parameter space, and is most accurate inside the "oceanographic funnel"
%  described in McDougall et al. (2003).  The GSW library function
%  "gsw_infunnel(SA,CT,p)" is avaialble to be used if one wants to test if
%  some of one's data lies outside this "funnel".
%
% INPUT:
%  SA   =  Absolute Salinity                                       [ g/kg ]
%  CT   =  Conservative Temperature (ITS-90)                      [ deg C ]
%  p    =  sea pressure                                            [ dbar ]
%           ( i.e. absolute pressure - 10.1325 dbar )
%  p_i  =  specific query points at which the interpolated SA_i and CT_i
%            are required                                          [ dbar ]
%
%  SA & CT need to have the same dimensions.
%  p may have dimensions Mx1 or 1xN or MxN, where SA & CT are MxN.
%  p_i needs to be either a vector or a matrix and have dimensions M_ix1
%  or M_ixN.
%
% OUTPUT:
%  SA_i = interpolated SA values at pressures p_i                  [ g/kg ]
%  CT_i = interpolated CT values at pressures p_i                 [ deg C ]
%
% AUTHOR:
%  Paul Barker and Trevor McDougall [ help@teos-10.org ]
%
% VERSION NUMBER: 3.06.1 (13th August, 2019)
%
% References
%  Barker, P.M., and T.J. McDougall, 2019: Two interpolation methods using
%   Multiply-Rotated Piecewise Cubic Hermite Interpolating Polynomials. 
%   J. Atmosph. Ocean. Tech. (Submitted).
%
%  McDougall, T.J., D.R. Jackett, D.G. Wright and R. Feistel, 2003:
%   Accurate and computationally efficient algorithms for potential
%   temperature and density of seawater.  J. Atmosph. Ocean. Tech., 20,
%   pp. 730-741.
%
%  Roquet, F., G. Madec, T.J. McDougall, P.M. Barker, 2015: Accurate
%   polynomial expressions for the density and specifc volume of seawater
%   using the TEOS-10 standard. Ocean Modelling., 90, pp. 29-43.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

[pl,number_of_profiles] = size(SA);
[mt,nt] = size(CT);
[mp,np] = size(p);
[interp_profile_length,np_i] = size(p_i);

if (pl ~= mt) | (number_of_profiles ~= nt)
    error('gsw_SA_CT_interp: SA and CT need to have the same dimensions')
end

if (mp == 1) & (np == 1)
    error('gsw_SA_CT_interp:  There must be at least 4 bottles')
elseif (number_of_profiles == np) & (mp == 1)
    p = p(ones(1,pl), :);
elseif (pl == mp) & (np == 1)
    p = p(:,ones(1,number_of_profiles));
elseif (number_of_profiles == mp) & (np == 1)
    p = p.';
    p = p(ones(1,pl), :);
elseif (pl == np) & (mp == 1)
    p = p.';
    p = p(:,ones(1,number_of_profiles));
elseif (pl == np) & (number_of_profiles == mp)
    p = p.';
elseif (pl == mp) & (number_of_profiles == np)
    % ok
else
    error('gsw_SA_CT_interp: Inputs array dimensions arguments do not agree')
end

if interp_profile_length == 1 & np_i > 1
    p_i = p_i.';
    dp_i = diff(p_i);
    if any(dp_i) < 0
        warning('gsw_SA_CT_interp: interpolating pressure must be monotonic')
        return
    end
    [interp_profile_length,np_i] = size(p_i);
elseif interp_profile_length == number_of_profiles & np_i~= number_of_profiles & all(diff(p_i,1,2)) >= 0
    p_i = p_i.';
    [interp_profile_length,np_i] = size(p_i);
elseif any(diff(p_i,1,1)) < 0
    warning('gsw_SA_CT_interp: interpolating pressure must be monotonic')
    return
else
    % Data shape and interval are ok.
end

%--------------------------------------------------------------------------
% Start of the calculation
%--------------------------------------------------------------------------

SA_i = NaN(interp_profile_length, number_of_profiles);
CT_i = SA_i;

for Iprofile = 1:number_of_profiles
    
    data_bottles = SA(:,Iprofile) + CT(:,Iprofile) + p(:,Iprofile);
    [Inn] = find(~isnan(data_bottles));
    if length(Inn) < 2
        SA_i(:,Iprofile) = NaN;
        CT_i(:,Iprofile) = NaN;
        continue
    end
    
    SA_obs = SA(Inn,Iprofile);
    CT_obs = CT(Inn,Iprofile);
    p_obs = p(Inn,Iprofile);
    dummy = 1e-3.*round(1e3.*p_obs);
    p_obs = dummy;
    if np_i > 1
        p_i_tmp = p_i(:,Iprofile);
    else
        p_i_tmp = p_i(:);
    end
    dummy = 1e-3.*round(1e3.*p_i_tmp);
    p_i_tmp = dummy;

    CTf_obs = gsw_CT_freezing_poly(SA_obs,p_obs);
    if any(CT_obs < (CTf_obs - 0.1))
        [Ifrozen] = find(CT_obs < (CTf_obs - 0.1));
        CT_obs(Ifrozen) = CTf_obs(Ifrozen);
    end
    
    pl = length(p_obs);
    Ishallow = 1:(pl-1);
    Ideep = 2:pl;
    dp_tmp = p_obs(Ideep) - p_obs(Ishallow);
    if any(dp_tmp <= 0)
        warning('gsw_SA_CT_interp: pressure must be monotonic')
        [p_sort,Ipsort] = (sort(p_obs));
        SA_sort = SA_obs(Ipsort);
        CT_sort = CT_obs(Ipsort);
        [p_obs, Ipunique] = unique(p_sort);
        SA_obs = SA_sort(Ipunique);
        CT_obs = CT_sort(Ipunique);
        pl = length(p_obs);
    end
    
    p_all = unique(sort([p_obs; p_i_tmp]));
    [Iobs_plus_interp] = find(p_all >= min(p_obs) & p_all <= max(p_obs));
    [Isurf_and_obs_plus_interp] = find(p_all <= max(p_obs));
    [dummy, Iout, I1] = intersect(p_i_tmp,p_all(Isurf_and_obs_plus_interp));
    [dummy, I2, I3] = intersect(p_obs,p_all(Iobs_plus_interp));
    
    clear interpolating_axis
    interpolating_axis(1:pl,1) = [0:pl-1];
    interpolating_axis_obs_plus_interp = gsw_pchip_interp(interpolating_axis,p_obs,p_all(Iobs_plus_interp));
    
    if ~exist('factor','var')
        factor = 9;
    end
    rec_factor = 1/factor;
    
    sin_pi_on_16 = 1.950903220161283e-1;  %sin_pi_on_16 = sin(pi./16)
    cos_pi_on_16 = 9.807852804032304e-1;  %cos_pi_on_16 = cos(pi./16)
    sin_pi_on_8 = 3.826834323650898e-1;   %sin_pi_on_8 = sin(pi./8)
    cos_pi_on_8 = 9.238795325112867e-1;   %cos_pi_on_8 = cos(pi./8)
    sin_3pi_on_16 = 5.555702330196022e-1; %sin_3pi_on_16 = sin(3pi./16)
    cos_3pi_on_16 = 8.314696123025452e-1; %cos_3pi_on_16 = cos(3pi./16)
    sin_pi_on_4 = 7.071067811865475e-1;   %sin_pi_on_8 = sin(pi./4)
    cos_pi_on_4 = 7.071067811865476e-1;   %cos_pi_on_8 = cos(pi./4)
    sin_5pi_on_16 = 8.314696123025452e-1; %sin_pi_on_8 = sin(5pi./16)
    cos_5pi_on_16 = 5.555702330196023e-1; %cos_pi_on_8 = cos(5pi./16)
    sin_3pi_on_8 = 9.238795325112867e-1;  %sin_pi_on_8 = sin(3pi./8)
    cos_3pi_on_8 = 3.826834323650898e-1;  %cos_pi_on_8 = cos(3pi./8)
    sin_7pi_on_16 = 9.807852804032304e-1; %sin_pi_on_8 = sin(7pi./16)
    cos_7pi_on_16 = 1.950903220161283e-1; %cos_pi_on_8 = cos(7pi./16)
   
    scaled_SA_obs = factor.*SA_obs;
    
    v1_tmp = CT_obs;
    q1_tmp = scaled_SA_obs;
    [v1_i,q1_i] = gsw_pchip_interp_SA_CT(v1_tmp,q1_tmp,interpolating_axis,interpolating_axis_obs_plus_interp);
    
    v2_tmp = scaled_SA_obs.*sin_pi_on_16 + CT_obs.*cos_pi_on_16;
    q2_tmp = scaled_SA_obs.*cos_pi_on_16 - CT_obs.*sin_pi_on_16;
    [v2_i,q2_i] = gsw_pchip_interp_SA_CT(v2_tmp,q2_tmp,interpolating_axis,interpolating_axis_obs_plus_interp);
        
    v3_tmp = scaled_SA_obs.*sin_pi_on_8 + CT_obs.*cos_pi_on_8;
    q3_tmp = scaled_SA_obs.*cos_pi_on_8 - CT_obs.*sin_pi_on_8;
    [v3_i,q3_i] = gsw_pchip_interp_SA_CT(v3_tmp,q3_tmp,interpolating_axis,interpolating_axis_obs_plus_interp);
        
    v4_tmp = scaled_SA_obs.*sin_3pi_on_16 + CT_obs.*cos_3pi_on_16;
    q4_tmp = scaled_SA_obs.*cos_3pi_on_16 - CT_obs.*sin_3pi_on_16;
    [v4_i,q4_i] = gsw_pchip_interp_SA_CT(v4_tmp,q4_tmp,interpolating_axis,interpolating_axis_obs_plus_interp);
    
    v5_tmp = scaled_SA_obs.*sin_pi_on_4 + CT_obs.*cos_pi_on_4;
    q5_tmp = scaled_SA_obs.*cos_pi_on_4 - CT_obs.*sin_pi_on_4;
    [v5_i,q5_i] = gsw_pchip_interp_SA_CT(v5_tmp,q5_tmp,interpolating_axis,interpolating_axis_obs_plus_interp);
        
    v6_tmp = scaled_SA_obs.*sin_5pi_on_16 + CT_obs.*cos_5pi_on_16;
    q6_tmp = scaled_SA_obs.*cos_5pi_on_16 - CT_obs.*sin_5pi_on_16;
    [v6_i,q6_i] = gsw_pchip_interp_SA_CT(v6_tmp,q6_tmp,interpolating_axis,interpolating_axis_obs_plus_interp);
        
    v7_tmp = scaled_SA_obs.*sin_3pi_on_8 + CT_obs.*cos_3pi_on_8;
    q7_tmp = scaled_SA_obs.*cos_3pi_on_8 - CT_obs.*sin_3pi_on_8;
    [v7_i,q7_i] = gsw_pchip_interp_SA_CT(v7_tmp,q7_tmp,interpolating_axis,interpolating_axis_obs_plus_interp);
      
    v8_tmp = scaled_SA_obs.*sin_7pi_on_16 + CT_obs.*cos_7pi_on_16;
    q8_tmp = scaled_SA_obs.*cos_7pi_on_16 - CT_obs.*sin_7pi_on_16;
    [v8_i,q8_i] = gsw_pchip_interp_SA_CT(v8_tmp,q8_tmp,interpolating_axis,interpolating_axis_obs_plus_interp);

    x = 0.125;
    
    % CT_i_1 = -q1_i.*sin(pi.*0/16) + v1_i.*cos(pi.*0/16);
    CT_i_1 = v1_i;
    % SA_i_1 = 0.25.*(q1_i.*cos(pi.*0/16) + v1_i.*sin(pi.*0/16));
    SA_i_1 = rec_factor.*q1_i;

    CT_i_2 = -q2_i.*sin_pi_on_16 + v2_i.*cos_pi_on_16;
    SA_i_2 = rec_factor.*(q2_i.*cos_pi_on_16 + v2_i.*sin_pi_on_16);
    
    CT_i_3 = -q3_i.*sin_pi_on_8 + v3_i.*cos_pi_on_8;
    SA_i_3 = rec_factor.*(q3_i.*cos_pi_on_8 + v3_i.*sin_pi_on_8);
    
    CT_i_4 = -q4_i.*sin_3pi_on_16 + v4_i.*cos_3pi_on_16;
    SA_i_4 = rec_factor.*(q4_i.*cos_3pi_on_16 + v4_i.*sin_3pi_on_16);
    
    CT_i_5 = -q5_i.*sin_pi_on_4 + v5_i.*cos_pi_on_4;
    SA_i_5 = rec_factor.*(q5_i.*cos_pi_on_4 + v5_i.*sin_pi_on_4);
    
    CT_i_6 = -q6_i.*sin_5pi_on_16 + v6_i.*cos_5pi_on_16;
    SA_i_6 = rec_factor.*(q6_i.*cos_5pi_on_16 + v6_i.*sin_5pi_on_16);
    
    CT_i_7 = -q7_i.*sin_3pi_on_8 + v7_i.*cos_3pi_on_8;
    SA_i_7 = rec_factor.*(q7_i.*cos_3pi_on_8 + v7_i.*sin_3pi_on_8);
    
    CT_i_8 = -q8_i.*sin_7pi_on_16 + v8_i.*cos_7pi_on_16;
    SA_i_8 = rec_factor.*(q8_i.*cos_7pi_on_16 + v8_i.*sin_7pi_on_16);
    
    CT_i_obs_plus_interp = x.*(CT_i_1 + CT_i_2 + CT_i_3 + CT_i_4 + CT_i_5 + CT_i_6 + CT_i_7 + CT_i_8);
    SA_i_obs_plus_interp = x.*(SA_i_1 + SA_i_2 + SA_i_3 + SA_i_4 + SA_i_5 + SA_i_6 + SA_i_7 + SA_i_8);
           
    [SA_i_limiting_obs_plus_interp, CT_i_limiting_obs_plus_interp] = gsw_linear_interp_SA_CT(SA_obs,CT_obs,interpolating_axis,interpolating_axis_obs_plus_interp);
        
    v_i_obs_plus_interp = gsw_specvol(SA_i_obs_plus_interp, CT_i_obs_plus_interp, p_all(Iobs_plus_interp));
    v_i_limiting = gsw_specvol(SA_i_limiting_obs_plus_interp, CT_i_limiting_obs_plus_interp, p_all(Iobs_plus_interp));
    
    [Ireplacenan] = find(isnan(v_i_obs_plus_interp) & ~isnan(v_i_limiting));
    SA_i_obs_plus_interp(Ireplacenan) = SA_i_limiting_obs_plus_interp(Ireplacenan);
    CT_i_obs_plus_interp(Ireplacenan) = CT_i_limiting_obs_plus_interp(Ireplacenan);
    v_i_obs_plus_interp(Ireplacenan) = v_i_limiting(Ireplacenan);
    
    SA_mid = 0.5*(SA_obs(1:end-1) + SA_obs(2:end));
    CT_mid = 0.5*(CT_obs(1:end-1) + CT_obs(2:end));
    p_mid = 0.5*(p_obs(1:end-1) + p_obs(2:end));
    
    v_shallower = gsw_specvol(SA_obs(1:end-1),CT_obs(1:end-1),p_mid);
    v_deeper = gsw_specvol(SA_obs(2:end),CT_obs(2:end),p_mid);
    v_mid = gsw_specvol(SA_mid,CT_mid,p_mid);
    delta_v_local = -v_mid + 0.5*(v_shallower + v_deeper);
    v_error = 2.*delta_v_local + 1e-7;
    v_error_obs_plus_interp = gsw_linear_interp(v_error,p_mid,p_all(Iobs_plus_interp));

    [max_v_obs, Imax_v_obs]  = max(v_i_obs_plus_interp(I3));
    max_v_data = max_v_obs + abs(v_error_obs_plus_interp(I3(Imax_v_obs)));
    if any(max(v_i_obs_plus_interp) > max_v_data)
        toolight = 1;
        while toolight == 1
            [Itoolight] = find(v_i_obs_plus_interp > max_v_data);
            [Ishallower] = find((I3 - Itoolight(1)) <= 0);
            Iabove = I2(Ishallower(end));
            Iabove_i = I3(Ishallower(end));
            if (Iabove+1) > I3(end)
                Ibelow_i = I3(end);
            else
                Ibelow_i = I3(Iabove + 1);
            end
            SA_i_obs_plus_interp(Iabove_i:Ibelow_i) = SA_i_limiting_obs_plus_interp(Iabove_i:Ibelow_i);
            CT_i_obs_plus_interp(Iabove_i:Ibelow_i) = CT_i_limiting_obs_plus_interp(Iabove_i:Ibelow_i);
            v_i_obs_plus_interp(Iabove_i:Ibelow_i) = v_i_limiting(Iabove_i:Ibelow_i);
            if ~any(max(v_i_obs_plus_interp) > max_v_data)
                toolight = 0;
            end
        end
    end
    
    [min_v_obs, Imin_v_obs]  = min(v_i_obs_plus_interp(I3));
    min_v_data = min_v_obs - abs(v_error_obs_plus_interp(I3(Imin_v_obs)));
    if any(min(v_i_obs_plus_interp) < min_v_data)
        tooheavy = 1;
        while tooheavy == 1
            [Itooheavy] = find(v_i_obs_plus_interp < min_v_data);
            [Ishallower] = find((I3 - Itooheavy(1)) <= 0);
            Iabove = I2(Ishallower(end));
            Iabove_i = I3(Ishallower(end));
            if (Iabove+1) > I3(end)
                Ibelow_i = I3(end);
            else
                Ibelow_i = I3(Iabove + 1);
            end
            SA_i_obs_plus_interp(Iabove_i:Ibelow_i) = SA_i_limiting_obs_plus_interp(Iabove_i:Ibelow_i);
            CT_i_obs_plus_interp(Iabove_i:Ibelow_i) = CT_i_limiting_obs_plus_interp(Iabove_i:Ibelow_i);
            v_i_obs_plus_interp(Iabove_i:Ibelow_i) = v_i_limiting(Iabove_i:Ibelow_i);
            if ~any(min(v_i_obs_plus_interp) < min_v_data)
                tooheavy = 0;
            end
        end
    end
    
    CTf_i_tointerp = gsw_CT_freezing_poly(SA_i_obs_plus_interp,p_all(Iobs_plus_interp));
    if any(CT_i_limiting_obs_plus_interp < (CTf_i_tointerp - 0.1))
        [ICTf_i_obs] = find(CT_i_limiting_obs_plus_interp < (CTf_i_tointerp - 0.1));
        CTf_i_tointerp(ICTf_i_obs) = CT_i_limiting_obs_plus_interp(ICTf_i_obs);
    end
    if any(CT_i_obs_plus_interp < (CTf_i_tointerp - 0.1))
        frozen = 1;
        while frozen == 1
            [Ifrozen] = find(CT_i_obs_plus_interp < (CTf_i_tointerp - 0.1));
            [Ishallower] = find((I3 - Ifrozen(1)) <= 0);
            Iabove = I2(Ishallower(end));
            Iabove_i = I3(Ishallower(end));
            if (Iabove+1) > I3(end)
                Ibelow_i = I3(end);
            else
                Ibelow_i = I3(Iabove + 1);
            end
            SA_i_obs_plus_interp(Iabove_i:Ibelow_i) = SA_i_limiting_obs_plus_interp(Iabove_i:Ibelow_i);
            CT_i_obs_plus_interp(Iabove_i:Ibelow_i) = CT_i_limiting_obs_plus_interp(Iabove_i:Ibelow_i);
            CTf_i_tointerp(Iabove_i:Ibelow_i) = gsw_CT_freezing_poly(SA_i_obs_plus_interp(Iabove_i:Ibelow_i),p_all(Iabove_i:Ibelow_i));
            if any(CT_i_limiting_obs_plus_interp < (CTf_i_tointerp - 0.1))
                [ICTf_i_obs] = find(CT_i_limiting_obs_plus_interp < (CTf_i_tointerp - 0.1));
                CTf_i_tointerp(ICTf_i_obs) = CT_i_limiting_obs_plus_interp(ICTf_i_obs);
            end
            if ~any(CT_i_obs_plus_interp < (CTf_i_tointerp - 0.1))
                frozen = 0;
            end
        end
    end
    
    [min_p_obs, Imin_p_obs] = min(p_obs);
    if min_p_obs ~= 0
        
        [Isurface] = find(p_i_tmp < min_p_obs);
        
        SA_i_tooutput = NaN(length(Isurf_and_obs_plus_interp),1);
        CT_i_tooutput = SA_i_tooutput;

        SA_i_tooutput(Isurface) = SA_i_obs_plus_interp(I3(Imin_p_obs));
        CT_i_tooutput(Isurface) = CT_i_obs_plus_interp(I3(Imin_p_obs));

        SA_i_tooutput(Iobs_plus_interp) = SA_i_obs_plus_interp;
        CT_i_tooutput(Iobs_plus_interp) = CT_i_obs_plus_interp;

    else
        
        SA_i_tooutput = SA_i_obs_plus_interp;
        CT_i_tooutput = CT_i_obs_plus_interp;

    end
    
    SA_i(Iout,Iprofile) = SA_i_tooutput(I1);
    CT_i(Iout,Iprofile) = CT_i_tooutput(I1);

end

end
