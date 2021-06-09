function[varargout]=maternimp(varargin)
%MATERNIMP  Impulse response function for the Matern random process.
%
%   [T,G]=MATERNIMP(N,ALPHA,H) returns first N points for the impulse 
%   response function G for a unit-variance complex-valued Matern process 
%   with slope parameter ALPHA and damping or range parameter H.
%
%   The impulse response function is also known as the Green's function.
%
%   T is an array of times T=[ESILON 1:1:N] where the first point is set to
%   a very small number, set to EPSILON=1e-6 to avoid the potential for 
%   infinite values at T=0 for some parameter settings.  
%
%   The sample interval is taken to be unity.  The damping parameter H is 
%   understood to have units of the inverse of the sample interval.
%
%   Note that for H=0, the case of fractional Brownian motion, G is   
%   formally defined but is not useful because it does not decay. 
%
%   The input parameters ALPHA and H may either be scalars or arrays of the 
%   same length M.  If the latter, then the output autocovariance function
%   R will be a matrix with N rows and M columns. 
%
%   See MATERNSPEC for a more thorough discussion of the Matern process.
%
%   For further details, see Sykulski, Olhede, Lilly, and Danioux (2015),
%   "Lagrangian time series models for ocean surface drifter trajectories."
%   __________________________________________________________________
%
%   Oscillatory Matern
%
%   MATERNIMP(N,ALPHA,C), where C is complex, uses H=REAL(C) for the
%   damping parameter and also sets a rotation frequency to OMEGA=-IMAG(C). 
%
%   The associated oscillatory process undergoes rotations at frequency 
%   OMEGA in addition to the damping at timescale 1/H.
%   __________________________________________________________________
%
%   Specified times
%
%   G=MATERNIMP(T,ALPHA,H) where the first argument is an *array*
%   rather than a scalar, will alternately return the impulse response at
%   the specified times T.  Note T is not output again in this case.
%   __________________________________________________________________
%
%   See also MATERNSPEC, MATERNCOV, MATERNOISE, MATERNPROPS, MATERNFIT.
%
%   'maternimp --t' runs some tests.
%
%   Usage:    [t,G]=maternimp(N,alpha,h);
%             G=maternimp(t,alpha,h);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2015  J.M. Lilly --- type 'help jlab_license' for details


 
if strcmpi(varargin{1}, '--t')
    maternimp_test,return
end

N=[];
tau=[];

if length(varargin)==3
    if length(varargin{1})==1
        N=varargin{1};
    else
        tau=varargin{1};
    end
end




if isempty(tau)
    tau=[0:1:N-1]'+1/2;
end

alpha=varargin{2};
H=real(varargin{3});
omegao=-imag(varargin{3});
arrayify(alpha,omegao,H);

G=zeros(size(tau,1),length(H));
if size(tau,2)~=length(H)
    tau=vrep(tau,length(H),2);
end
for i=1:length(H)
    U=frac(1,2).*(1+sign(tau(:,i)));
    d=sqrt(frac(2*H(i).^(2*alpha(i)-1),maternc(alpha(i))));
    G(:,i)=d.*U.*(tau(:,i).^(alpha(i)-1)).*exp(-tau(:,i).*H(i)).*exp(sqrt(-1)*omegao(i)*tau(:,i)).*frac(1,gamma(alpha(i)));
end

if ~isempty(N)
    varargout{1}=tau;
    varargout{2}=G;
else
    varargout{1}=G;
end



function[]=maternimp_test
alpha=[1 1.5 2 3 4];
h=[.01 .02 .05 .2 1]/10;
[alpha,h]=meshgrid(alpha,h);
vcolon(alpha,h);
[tau,g]=maternimp(15000,alpha,h);
d=sqrt(frac(2*h.^(2*alpha-1),maternc(alpha)));
for i=1:size(g,2)
    g(:,i)=g(:,i)./d(i);
end
reporttest('MATERNIMP integral of Green''s function matches analytic form',aresame((sum(g,1)'.*(h.^alpha))'-1,zeros(size(h')),2.5e-2))

t=[0:1:150]';
alpha=1.4;h=1/10;
g=maternimp(t,alpha,h);
R1=conv(g,flipud(g));
[tau,R2]=materncov(length(t),1,alpha,h,'full');
eps=sum(squared(R1-R2))./sum(squared(R1));
reporttest('MATERNIMP Green''s function convolves to autocovariance sample rate',eps<1e-3)


N=1000;
alpha=1.4;h=1/10;
[tau,g]=maternimp(N,alpha,h);
R1=conv(g,flipud(g));
[tau,R2]=materncov(N,1,alpha,h,'full');
eps=sum(squared(R1-R2))./sum(squared(R1));
reporttest('MATERNIMP Green''s function convolves to autocovariance without time input',eps<1e-3)

