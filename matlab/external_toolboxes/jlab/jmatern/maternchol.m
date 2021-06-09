function[T]=maternchol(N,alpha,H,alg)
%MATERNCHOL  Cholesky decomposition of the Matern covariance matrix.
%
%   MATERNCHOL is a low-level function called by MATERNOISE.
%
%   T=MATERNCHOL(N,ALPHA,H) returns the Cholesky decomposition of the N x N
%   autocovariance matrix of a unit variance complex-valued Matern process
%   with slope parameter ALPHA and damping or range parameter H.  
%
%   Note that H=0 case also works.  This corresponds to fractional Brownian
%   motion, whereas a Matern process is formally defined for H > 0.  In
%   this case, variance is undefined, and autovariance matrix is for a 
%   unit spectral amplitude process rather than a unit variance process.
%  
%   See MATERNSPEC for a more thorough discussion of the Matern process.
%  
%   For further details, see Sykulski, Olhede, Lilly, and Danioux (2015),
%   "Lagrangian time series models for ocean surface drifter trajectories."
%   __________________________________________________________________
% 
%   Oscillatory Matern
%
%   MATERNCHOL(N,ALPHA,C), where C is complex, uses H=REAL(C) for the
%   damping parameter and also sets a rotation frequency to OMEGA=-IMAG(C). 
%  
%   The associated oscillatory process undergoes rotations at frequency 
%   OMEGA in addition to the damping at timescale 1/H.
%   __________________________________________________________________
%
%   'maternchol --t' runs a test.
%
%   Usage: T=maternchol(N,alpha,h);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2015  A.M. Sykulski and J.M. Lilly
%                                 --- type 'help jlab_license' for details 
if strcmp(N, '--t')
    maternchol_test,return
end

if nargin == 3
    alg='fast';
end

omegao=-imag(H);
H=real(H);

if H==0
    C=fbm_matrix(N,alpha,alg);
else
    C=matern_matrix(N,alpha,omegao,H,alg);
end

T=zeros(size(C));
%Note lower triangular is a causal filter.  Plot rows (not columns).
%T=chol(conj(C(:,:,k)),'lower');
%This is because of issues with non-positive definiteness, see chol
[temp,p]=chol(C,'lower');
T(1:size(temp,1),1:size(temp,2))=conj(temp);


function [matrix]=matern_matrix(N,alpha,omegao,h,alg)

M=length(alpha);

if strcmpi(alg(1:3),'loo')
    %This is from Adam's original code, used only for testing purposes;
    %does not work when omegao is nonzero or for non-unit sample rate
    tpz = zeros(N,1);
    d=frac(2*pi *H.^(2*alpha-1),beta(1/2,alpha-1/2));
    tpz(1)=d.*((sqrt(pi))/((2^(alpha-1.5))*gamma(alpha)*(h^(2*alpha-1))))*gamma(alpha-1/2)*pow2(alpha-3/2);
    for n=1:N-1
        tpz(n+1)=d.*((sqrt(pi))/((2^(alpha-1.5))*gamma(alpha)*(h^(2*alpha-1))))*((h*n)^(alpha-.5))*besselk(alpha-.5,h*n);
    end
    tpz=tpz./(2*pi);
    matrix = toeplitz(tpz);
else
    %My translation
    [tau,tpz]=materncov(N,1,alpha,h-sqrt(-1)*omegao);
    %matrix=jtoeplitz(tpz);
    %figure,plot(tpz)
    matrix=toeplitz(tpz);
end


function [matrix]=fbm_matrix(N,alpha,alg)

% This is coded up as in Barton and Poor (1988) 
%   "Signal Detection in Fractional Gaussian Noise"
% IEEE TRANSACTIONS ON INFORMATION THEORY, VOL. 34, NO. 5, SEPTEMBER 1988

% Note the Hurst parameter H = alpha  - 1/2
alpha(alpha==0.5)=0.5+1e-5;
alpha(alpha==1.5)=1.5-1e-5;
alpha(alpha==1)=1+1e-5;

%H = alpha - 1/2-sqrt(-1)*omegao;
H = alpha - 1/2;
M=length(alpha);
t=[1:N]';

if strcmpi(alg(1:3),'loo')
    matrix=zeros(N,N);
    
    %This is from Adam's original code, used only for testing purposes;
    %does not work when omegao is nonzero
    for j=1:N
        for k=1:N
            matrix(j,k)=(1./2)...
                .*(abs(t(j)).^(2*H)+abs(t(k)).^(2*H)-abs(t(k)-t(j)).^(2*H));
        end
    end
    matrix=matrix.*(-gamma(2-2*H).*cos(pi*H))./(pi.*H.*(2*H-1));
else
    %My translation
    fact=frac(1,2).*frac(-gamma(2-2*real(H)).*cos(pi.*H),pi.*H.*(2.*H-1));
    tk2H=vrep(realpow(t,2*H),N,2);
    tk=vrep(t,size(t,1),2);
    dtk2H=realpow(abs(tk-tk'),2*H);
    matrix=fact*(tk2H+tk2H'-dtk2H);
end




function[]=maternchol_test


N = 100; 
alpha=1.6;
h=0.15;

t1=tic;C1 = matern_matrix(N,alpha,0,h,'chol');etime1=toc(t1);
t2=tic;C2 = matern_matrix(N,alpha,0,h,'fast');etime2=toc(t2);

reporttest('MATERNCHOL loop and loopless algorithms match for Matern', aresame(C1./C2,ones(size(C1)),1e-10))
%disp(['Loopless algorithm was ' num2str(etime1./etime2) ' times faster than loop.'])

t1=tic;C1 = fbm_matrix(N,alpha,'chol');etime1=toc(t1);
t2=tic;C2 = fbm_matrix(N,alpha,'fast');etime2=toc(t2);

reporttest('MATERNCHOL loop and loopless algorithms match for fBM', aresame(C1./C2,ones(size(C1)),1e-10))
%disp(['Loopless algorithm was ' num2str(etime1./etime2) ' times faster than loop.'])

%alpha=0.505+rand(1,N/2)/1.1;%min(alpha),max(alpha)
alpha=1+rand(1,N/2)/1.1;%min(alpha),max(alpha)
h=rand(1,N/2);

[C1,C2]=vzeros(N,N,length(alpha));
t1=tic;
for i=1:length(alpha)
    C1(:,:,i) = matern_matrix(N,alpha(i),0,h(i),'chol');
end
etime1=toc(t1);
   
t2=tic;
for i=1:length(alpha)
    C2(:,:,i) = matern_matrix(N,alpha(i),0,h(i),'fast');
end
etime2=toc(t2);

reporttest('MATERCHOL loop and loopless algorithms match for Matern with array input', aresame(C1./C2,ones(size(C1)),1e-10))
%disp(['Loopless algorithm was ' num2str(etime1./etime2) ' times faster than loop.'])

alpha(alpha<=0.5)=0.5+1e-5;
alpha(alpha>=1.5)=1.5-1e-5;

[C1,C2]=vzeros(N,N,length(alpha));
t1=tic;
for i=1:length(alpha)
    C1(:,:,i) =  fbm_matrix(N,alpha(i),'chol');
end
etime1=toc(t1);
t2=tic;
for i=1:length(alpha)
    C2(:,:,i) = fbm_matrix(N,alpha(i),'fast');
end
etime2=toc(t2);

reporttest('MATERCHOL loop and loopless algorithms match for fBm with array input', aresame(C1./C2,ones(size(C1)),1e-10))
%disp(['Loopless algorithm was ' num2str(etime1./etime2) ' times faster than loop.'])

%T=maternchol(1000,1.5,1/100);

N = 500;
alpha1=[1:.1:3.5];
h1=logspace(-1,0,20)/5;
[alpha,h]=meshgrid(alpha1,h1);

[G,G1,G2]=vzeros(length(h1),length(alpha1),N);
T=zeros(length(h1),length(alpha1),N,N);
for j=1:length(alpha1)
    for i=1:length(h1);
        Ttemp=maternchol(N,alpha1(j),h1(i));
        G(i,j,:)=flipud(Ttemp(end,:)');
        [t,G1(i,j,:)]=maternimp(N,alpha1(j),h1(i));  
    end
end

err=frac(sum(squared(G-G1),3),sum(squared(G),3));
reporttest('MATERCHOL matches form of impulse response function to within 1%', allall(err<0.01))
%See MATERNOISE for computation of error term.
