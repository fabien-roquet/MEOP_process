function [Tcorr,Scorr,Thp] = thermal_cell_correction (alpha,beta,dt,S,T,P,direction)

% input:
% alpha: =[alpha_T alpha_C], initial error for T and C, size([2])
% beta: inverse relaxation pressure interval dbar-1, size([1])
% dt: sampling interval in dbar, size([1])
% P: Pressure from tag (interpolated on a 1dbar vertical grid, dbar), size [nr np]
% T: Temperature from tag (interpolated on a 1dbar vertical grid, degC), size [nr np]
% S: Salinity from tag (interpolated on a 1dbar vertical grid, psu), size [nr np]
% direction [default: 1]: = 1 if upcast, 0 if downcast
%
% output:
% Tcorr: corrected T
% Scorr: corrected S
% Thp  : hig-passed temperature signal

[nr np] = size(S);
if nr < 3, Scorr = S; return; end
if ~exist('direction','var'), direction=1; end

C = gsw_C_from_SP(S(:),T(:),P(:));
C = reshape(C,size(S));

if direction
    P1 = flipud(P);
    T1 = flipud(T);
    S1 = flipud(S);
    C1 = flipud(C);
else
    P1 = P;
    T1 = T;
    S1 = S;
    C1 = C;
end

tau = 1/beta;
R = tau/(tau+dt);
Thp = zeros(size(S1));
Thp(1,:)=0;
for kk=2:nr;
    Thp(kk,:) = R * ( Thp(kk-1,:) + T1(kk,:) - T1(kk-1,:) );
    Thp(kk,isnan(T1(kk-1,:))) = 0;
end

gradC = gsw_C_from_SP(S1(:),T1(:)+0.5,P1(:))-gsw_C_from_SP(S1(:),T1(:)-0.5,P1(:));
gradC = reshape(gradC,size(P1));
gradC(isnan(gradC)) = 0;

Tc = alpha(1)/(1-.5*beta*dt)       .*Thp;
Cc = alpha(2)/(1-.5*beta*dt).*gradC.*Thp;

Tcorr = T1 + Tc;
Ccorr = C1 + Cc;
Scorr = gsw_SP_from_C(Ccorr(:),Tcorr(:),P1(:));
Scorr = reshape(Scorr,size(P1));

if direction
    Thp   = flipud(Thp);
    Tcorr = flipud(Tcorr);
    Ccorr = flipud(Ccorr);
    Scorr = flipud(Scorr);
end

