function[varargout]=materncov(varargin)
%MATERNCOV  Autocovariance function of the Matern random process.
%
%   [TAU,R]=MATERNCOV(N,SIG,ALPHA,H) gives the theoretical autocovariance 
%   function for a length N complex-valued Matern process with variance
%   SIG^2, slope parameter ALPHA, and damping or range parameter H.
%
%   TAU is an array of time lags at which R is computed, and is given by 
%   TAU=[0,1,...,N-1].
%
%   By definition, R is one-sided theoretical autocovariance at 
%   non-negative time lags.  See below for the relationship between this 
%   and the full, length (2N-1) theoretical autocovariance function. 
%
%   The sample interval is taken to be unity.  The damping parameter H is 
%   understood to have units of the inverse of the sample interval.
%
%   Note that for H=0, the case of fractional Brownian motion, R will 
%   contain only INFs because the autocovariance function is unbounded.
%
%   The input parameters SIG, ALPHA, and H, may all either be scalars or 
%   arrays of the same length M.  If the latter, then the output 
%   autocovariance function R will be a matrix with N rows and M columns. 
%
%   See MATERNSPEC for a more thorough discussion of the Matern process.
%
%   For further details, see Sykulski, Olhede, Lilly, and Danioux (2015),
%   "Lagrangian time series models for ocean surface drifter trajectories."
%   __________________________________________________________________
%
%   Relationship to full autocovariance
%
%   For a time series of length N, the full autocovariance function RF is 
%   length 2N-1, defined at time lags -N+1,-N+2...,-1,0,1,...,N-2,N-1.
%
%   The one-sided autocovariance R contains the full autocovariance RF at 
%   positive time lags. Negative lags are given by Hermitian symmetry.
%
%   [TAUF,RF]=MATERNCOV(...,'full') returns the full (two-sided)
%   autocovariance RF and the corresponding two-sided time array TAUF. 
%
%   RF is constructed from R as RF=[FLIPUD(CONJ(R(2:end,:));R]. 
%   __________________________________________________________________
%
%   Oscillatory Matern
%
%   MATERNCOV(N,SIG,ALPHA,C), where C is complex, uses H=REAL(C) for the
%   damping parameter and also sets a rotation frequency to OMEGA=-IMAG(C). 
%
%   The associated oscillatory process undergoes rotations at frequency 
%   OMEGA in addition to the damping at timescale 1/H.
%   __________________________________________________________________
%
%   Specified times
%
%   R=MATERNCOV(TAU,SIG,ALPHA,H) where the first argument is an *array*
%   rather than a scalar, will alternately return the autocovariance at 
%   the specified time lags TAU.  TAU is not output again in this case.
%   __________________________________________________________________
%
%   See also MATERNSPEC, MATERNIMP, MATERNOISE, MATERNPROPS, MATERNFIT.
%
%   'materncov --t' runs some tests.
%
%   Usage:    [tau,R]=materncov(N,sig,alpha,h);
%             R=materncov(tau,sig,alpha,h);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2015  J.M. Lilly --- type 'help jlab_license' for details

if strcmpi(varargin{1}, '--t')
    materncov_test,return
end

sid='one';

N=[];
tau=[];

if ischar(varargin{end})
    sid=varargin{end};
    varargin=varargin(1:end-1);
end

btauin=false;

if length(varargin{1})==1
    N=varargin{1};
else
    btauin=true;
    tau=varargin{1};
end

if isempty(tau)
    tau=[0:N-1]';
end

sig=varargin{2};
alpha=varargin{3};
H=real(varargin{4});
omegao=-imag(varargin{4});
arrayify(sig,alpha,omegao,H);
R=zeros(size(tau,1),length(H));

for i=1:length(H)
    if sig(i)==0
        R(:,i)=0;
    else
        R(:,i)=sig(i).^2.*maternfun(alpha(i)-1/2,H(i).*abs(tau));
        if omegao(i)~=0
            R(:,i)=R(:,i).*exp(sqrt(-1)*abs(tau)*omegao(i));
        end
    end
end

if btauin
    varargout{1}=R;
else
    if strcmpi(sid(1:3),'ful')
        tau=[-flipud(tau(2:end,:));tau];
        R=[flipud(conj(R(2:end,:)));R];
    end
    varargout{1}=tau;
    varargout{2}=R;
end

function[]=materncov_test
alpha=1.5;
h=0.1;
sig=1;

N=1000;
[tau,R]=materncov(N,sig,alpha,h);
[f,Spp,Snn]=maternspec(N,sig,alpha,h);
R(1)=R(1)/2;
S=2*real(fft(R));
Spp2=S(1:length(Spp));
S=flipud(S);
Snn2=[S(end);S(1:length(Snn)-1)];

b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-4);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),1e-4);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC for case with negligible aliasing, even N',b1&&b2)

N=1000-1;
[tau,R]=materncov(N,sig,alpha,h);
R(1)=R(1)/2;
[f,Spp,Snn]=maternspec(N,sig,alpha,h);
S=2*real(fft(R));
Spp2=S(1:length(Spp));
S=flipud(S);
Snn2=[S(end);S(1:length(Snn)-1)];

b1=aresame(Spp2./maxmax(Spp),Spp./maxmax(Spp),1e-4);
b2=aresame(Snn2./maxmax(Snn),Snn./maxmax(Snn),1e-4);

reporttest('MATERNCOV Fourier transforms to MATERNSPEC for case with negligible aliasing, frequency shift case, odd N',b1&&b2)


% alpha=1.5;
% h=0.1;
% sig=1;
% omegao=0.1;
% 
% N=1000;
% figure
% [tau,R]=materncov(N,sig,1,h-1i*omegao);
% [f,Spp,Snn]=blurspec(R);
% plot(-f,Snn),hold on,plot(f,Spp)
% [tau,R]=materncov(N,sig,1,h+1i*omegao);
% [f,Spp,Snn]=blurspec(R);
% plot(-f,Snn),hold on,plot(f,Spp)
%  
% [f,Spp,Snn]=maternspec(N,sig,alpha,h);

% N=1000;
% r1=[0:.00001:0.00001]';
% alpha=[1/2:0.1:4];
% R=materncov(r1,1,alpha(:),1/100);
% dr=(R(2,:)-R(1,:))./(r1(2)-r1(1));
% dr2=-frac(gamma(abs(alpha-3/2)),gamma(alpha-1/2)).*frac(1/100,2);
% %dr2(alpha==1)=-1/100;
% 
% plot(alpha,-dr,'r*'),hold on,plot(alpha,-dr2,'o')
% 
% dxi=-frac(r,2).*frac(gamma(alpha-1),gamma(alpha)).*maternfun(alpha-1,r);
% reporttest('MATERNFUN iterative expression for derivative matches direct expression',aresame(xi1(2:end,:),dxi(2:end,:),1e-10))
% 

