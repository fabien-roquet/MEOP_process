function[xi,xi1,xi2]=maternfun(alpha,r)
%MATERNFUN  Returns the Matern function.
%
%   MATERNFUN is a low-level function called by MATERNCOV.
%
%   XI=MATERNFUN(ALPHA,R) returns the value of the Matern function of order
%   ALPHA at location R.  THe Matern function is defined as
%
%      XI(R)=1/(GAMMA(ALPHA)*2^(ALPHA-1)) * R^ALPHA K_ALPHA(R)
%   
%   where K_ALPHA is the modified Bessel function of the second kind of 
%   order ALPHA. This definition implies XI(0)=1.
%
%   ALPHA and R can either be scalars or 1D arrays.  XI will be an array of
%   size LENGTH(R) x LENGTH(ALPHA).
%
%   [XI,XI1,XI2]=MATERNFUN(ALPHA,R) also returns its first two derivatives.
%
%   'maternfun --t' runs some tests.
%   'maternfun --f' generates a sample figure.
%
%   See also MATERNRAD.
%
%   Usage:  xi=maternfun(alpha,r);
%          [xi,xi1,xi2]=maternfun(alpha,r);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2015 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmpi(alpha, '--t')
    maternfun_test,return
end

% if size(alpha,1)>1
%     error('MATERNFUN is expecting a scalar or a row vector for ALPHA.')
% end

r=abs(r);
%r(abs(r)<1e-10)=1e-10;
%[alpha,r]=arrayify(alpha(:)',r(:));

fact=frac(1,gamma(alpha).*pow2(alpha-1));

if length(alpha)==1
    if alpha~=0
        xi=fact.*(r.^alpha).*besselk(alpha,r);
    else
        xi=fact.*exp(-r);
    end
else
    xi=fact.*(r.^alpha).*besselk(alpha,r);
    xi(alpha==0)=(fact(alpha==0).*exp(-r(alpha==0)));
end
    
xi(r==0)=1;   %fix for numerical problem at r=0
%See 9.6.9 of Abramowitz and Stegun for the small-argument behavior of K


%xi=real(xi);
if nargout>1
    xi1=-fact.*(r.^alpha).*besselk(alpha-1,r);
end
if nargout>2
    xi2=fact.*(r.^alpha).*besselk(alpha-2,r)-fact.*(r.^(alpha-1)).*besselk(alpha-1,r);
end
 
function[]=maternfun_test
alpha=[1/2:1/4:10];
xio=maternfun(alpha,1e-10);
reporttest('MATERNFUN anticipated value at x=0',aresame(xio,1+0*xio,1e-8))

%alpha=[1/2+1/10 3/4:1/4:10];
%[xi,xi1,xi2]=maternfun(alpha,maternrad(alpha));
%reporttest('MATERNFUN second derivative vanishes at MATERNRAD',xi2(1)./xi(1)<1e-2&&allall(xi2(2:end)./xi(2:end)<1e-3))

r1=[0:.01:10]';
alpha1=[1/2:0.1:4];
[alpha,r]=meshgrid(alpha1,r1);
[xi,xi1]=maternfun(alpha,r);
dxi=vdiff(xi,1)./0.01;
reporttest('MATERNFUN iterative expression for derivative matches derivative',aresame(xi1(2:end,:),dxi(2:end,:),0.03))

r1=[0:.01:10]';
alpha1=[1/2:0.1:0.9 1.1:0.1:4];
[alpha,r]=meshgrid(alpha1,r1);
[xi,xi1]=maternfun(alpha,r);
dxi=-frac(r,2).*frac(gamma(alpha-1),gamma(alpha)).*maternfun(alpha-1,r);
reporttest('MATERNFUN iterative expression for derivative matches direct expression',aresame(xi1(2:end,:),dxi(2:end,:),1e-10))

function[]=maternfun_fig
r1=[0:.01:10]';
alpha1=[1/2:0.1:4];
[alpha,r]=meshgrid(alpha1,r1);
xi=maternfun(alpha,r);
figure,plot(r1,xi)
hold on,h=plot(r1,xi(:,1)); linestyle -h h 3D
hold on,h=plot(r1,xi(:,end)); linestyle -h h 3k
title('The Matern function from \alpha = 1/2 (gray) to \alpha = 4 (black)')




