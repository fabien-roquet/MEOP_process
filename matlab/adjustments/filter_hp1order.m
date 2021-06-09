function Xhp = filter_hp1order (tau,dt,X)

% input:
% tau: relaxation time, size([1])
% dt: sampling interval in s, size([1])
% X: input signal, size [nr np], nr number of samples, np number of profiles
% Xhp: high-passed filtered signal

[nr np] = size(X);
if nr < 3, Xhp = X*0; return; end

R = tau/(tau+dt);
Xhp = X*0;
Xhp(1,:) = 0;
for kk=2:nr;
    Xhp(kk,:) = R * ( Xhp(kk-1,:) + X(kk,:) - X(kk-1,:) );
    Xhp(kk,isnan(X(kk-1,:))) = 0;
end

