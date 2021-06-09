function Scorr = thermal_cell_correction (alpha,beta,dt,S,T,P)

% input:
% alpha: initial error, size([1])
% beta: inverse relaxation time, size([1])
% dt: sampling interval in dbar, size([1])
% P: Pressure from tag (interpolated on a 1dbar vertical grid), size [nr np]
% T: Temperature from tag (interpolated on a 1dbar vertical grid), size [nr np]
% S: Salinity from tag (interpolated on a 1dbar vertical grid), size [nr np]
%
% Scorr: corrected S

[nr np] = size(S);
if nr < 3, Scorr = S; return; end
if ~exist('direction','var'), direction=1; end

C = gsw_C_from_SP(S(:),T(:),P(:));
C = reshape(C,size(S));

P1 = flipud(P);
T1 = flipud(T);
S1 = flipud(S);
C1 = flipud(C);

tau = 1/beta;
R = tau/(tau+dt);
Thp = zeros(size(S1));
Thp(1,:)=0;
for kk=2:nr;
    Thp(kk,:) = R * ( Thp(kk-1,:) + T1(kk,:) - T1(kk-1,:) );
    Thp(kk,isnan(T1(kk-1,:))) = 0;
end

gradS = gsw_SP_from_C(C1,T1+0.5,P1)-gsw_SP_from_C(C1,T1-0.5,P1);
gradS(isnan(gradS)) = 0;

Sc = alpha/(1-.5*beta*dt).*gradS.*Thp;
Scorr = S1 + Sc;

Scorr = flipud(Scorr);

