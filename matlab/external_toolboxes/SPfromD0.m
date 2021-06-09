function SA = SPfromD0(D0,PT)
% polynomial fit of the inverse of sw_dens0

x_SPfromD0 = [
    -16964.6978236309
    -30713.2365330974
    8831.49422771523
    -30.1672464989239
    0.0317617735630268
    -1.03949918913186e-05
    568.509407564307
    -1225.72005026648
    -60.7199732371302
    0.102240530716968
    -4.04385345442789e-05
    1081.07839356717
    16.9734428585712
    -0.00317527946514473
    3.42199986752477e-06
    1.27970501600622
    -0.258026403569788
    0.000238030779684243
    0.812050370313116
    4.32426768533235e-05
    -0.00686852290612419
    -16.9595580646074
    -29.4566335010399
    1.81133526466634
    0.793421792097501
    9.95159324377875
    60.7269655332568
    ];
N1=5;
N2=2;

Xs1=zeros(N1+1,N1+1,1); nn1=(N1+1)*(N1+2)/2;
Xs2=zeros(N2+1,N2+1,1); nn2=(N2+1)*(N2+2)/2;

Ioptim=zeros(N1+1,N1+1,1); nn1=(N1+1)*(N1+2)/2;
for kk=1:N1+1,
    Ioptim(kk,1:N1+2-kk,1)=1;
end
Ioptim1=find(Ioptim);

Ioptim=zeros(N2+1,N2+1,1);
for kk=1:N2+1,
    Ioptim(kk,1:N2+2-kk,1)=1;
end
Ioptim2=find(Ioptim);
Ioptim3=[Ioptim1;Ioptim2];

Xs1(Ioptim3(1:nn1))=x_SPfromD0(1:nn1);
Xs2(Ioptim3(nn1+(1:nn2)))=x_SPfromD0(nn1+(1:nn2));

SA = compute_poly(D0,PT,D0*0,Xs1)./compute_poly(D0,PT,D0*0,Xs2);




%%
function v = compute_poly(S,T,Z,X)
% function [s,a,b]=compute_poly(S,T,Z,X,s0,t0,z0)
%  compute polynomial values 
%  input data: S, T, Z: sali, temp and depth
%              X=[X]_{klm}: polynomial coeff (3-dim array)
%
% v = sum( X(k,l,m) * S^k * T^l * Z^m )
v = S*0;
for kk=size(X,1)-1:-1:0,
    aux = v*0;
    for ll=size(X,2)-1:-1:0,
        p=flipud(squeeze(X(kk+1,ll+1,:)));
        if all(p==0),
            aux = aux .* T ;
        else
            while p(1)==0, p=p(2:end); end
            aux = aux .* T + polyval(p,Z);
        end
    end
    v = v .* S + aux ;
end



