function [Salcor,Tc]=sal_cor (alpha,beta,SI,T,S,P)

% input:
% alpha: initial error size([1 1])
% beta: inverse relaxation time size([1 1])
% SI: sampling interval in s: size([1 1])
% T: Temperature from tag size ([m 1]) or ([1 m])
% S: Salinity from tag , same size as T
% P: Pressure from tag, same size as T



Ptemp=P;
if Ptemp(1)<Ptemp(end)
    P=P(end:-1:1);
    T=T(end:-1:1);
    S=S(end:-1:1);
end


dt=T(2:end)-T(1:end-1); dt(isnan(dt))=0;
T68=1.00024.*T;
Ctag=sw_cndr(S,T68,P).*sw_c3515;
Sal_plus=sw_salt(Ctag/sw_c3515,T68+0.5,P);
Sal_min=sw_salt(Ctag/sw_c3515,T68-0.5,P);
gradS=Sal_plus-Sal_min; gradS(isnan(gradS))=-1;

a=(2*alpha)./((SI*beta)+2);
b=1-((2*a)/alpha);
Salcor=zeros(length(S),1);
Scini=0;
Sc(2)=-b*Scini+a*(gradS(2)*dt(1));%

Tcini=0;
Tc(2)=-b*Tcini+dt(1);%

Salcor(1)=S(1);
Salcor(2)=S(2)+Sc(2);

for k=3:length(S);
    if P(k)-P(k-1)>3
        Tc(k)=0;
        Sc(k)=0;
        Salcor(k)=S(k);
        dt(k-1)=dt(k);
    else
        Tc(k)=-1.*b.*Tc(k-1)+dt(k-1);
        Sc(k)=-1.*b.*Sc(k-1)+a*(gradS(k).*dt(k-1));
        Salcor(k)=S(k)+Sc(k);
    end
end

if Ptemp(1)<Ptemp(end)
    Tc=Tc(end:-1:1);
    Salcor=Salcor(end:-1:1);
end

