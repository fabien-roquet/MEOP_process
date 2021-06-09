function[varargout]=maternspec(varargin)
%MATERNSPEC  Fourier spectrum of the Matern random process.
% 
%   [F,S]=MATERNSPEC(N,SIG,ALPHA,H) returns the spectrum S of a length N
%   complex-valued Matern random process having variance SIG^2, slope 
%   parameter ALPHA, and range or damping parameter H.
%
%   F is an array of one-sided (positive) Fourier frequencies for a time
%   series of length N, F=FOURIER(N).  Note that F is a *radian* frequency. 
%   For simplicity, the sample interval is taken to be unity.
%
%   The lengths of the output variables F and S are N/2+1 for even N, and
%   (N+1)/2 for odd N.
%     
%   S is the postive rotary spectrum given by
%
%        S(F) = SIG^2 / (F^2 + H^2)^ALPHA * D(ALPHA,H)
%    
%   where D(ALPHA,H)= 2 * H^(2*ALPHA-1)./MATERNC(ALPHA) is a normalizing
%   function, with MATERNC being a function of the slope parameter ALPHA. 
%
%   Note that the positive and negative rotary spectra are identical for
%   the standard Matern process.
%
%   The input parameters SIG, ALPHA, and H, may all either be scalars or 
%   arrays of the same length M.  If the latter, then the output spectrum
%   S will be a matrix with LENGTH(F) rows and M columns. 
%
%   The sample interval is taken to be unity.  The damping parameter H is 
%   understood to have units of the inverse of the sample interval.
%
%   For further details, see Sykulski, Olhede, Lilly, and Danioux (2015),
%   "Lagrangian time series models for ocean surface drifter trajectories."
%   __________________________________________________________________
%
%   Sample interval  
%
%   MATERNSPEC(DT,N,SIG,ALPHA,H) uses DT as the sample interval, instead of
%   the default DT=1.  H is understood to have the same units as 1/DT.
%   __________________________________________________________________
%
%   Oscillatory Matern
%
%   [F,SPP,SNN]=MATERNSPEC(N,SIG,ALPHA,C), with C complex, uses H=REAL(C)
%   for the damping, and also sets a rotation frequency to OMEGA=-IMAG(C).  
%
%   SPP and SNN are now the postive rotary and negative rotary spectra,
%
%        SPP(F) = SIG^2 / [(F-OMEGA)^2 + H^2]^ALPHA    * D(ALPHA,H)
%        SPP(F) = SIG^2 / [(-F-OMEGA)^2 + H^2]^ALPHA   * D(ALPHA,H)
%
%   with the spectrum for positive frequencies +F returned in SPP, and for 
%   negative frequencies -F in SNN.  
%
%   The associated oscillatory process undergoes rotations at frequency 
%   OMEGA in addition to the damping at timescale 1/H.  
%
%   Note that OMEGA has units of radians per sample interval.  If DT is 
%   input, then OMEGA is understood to ahve the same units as 1/DT.
%   __________________________________________________________________
%   
%   Special cases
%
%   The generalized Matern family contains several other families of
%   random processes as special cases:
%
%     OMEGA=0           - The complex-valued Matern process
%     OMEGA=0 and H=0   - Complex-valued fractional Brownian motion
%     ALPHA=1           - The complex Ornstein-Uhlenbeck process
%     ALPHA=1, OMEGA=0  - Z=X+iY where X and Y are real-valued OU processes 
%
%   See Sykulski, Olhede, Lilly, and Danioux (2015) for more discussion.
%   __________________________________________________________________
%
%   See also MATERNCOV, MATERNOISE, MATERNPROPS, MATERNFIT, BLURSPEC.
%
%   'maternspec --f' generates some sample figures.
%
%   Usage:      [f,s]=maternspec(N,sig,alpha,h);
%               [f,spp,snn]=maternspec(dt,N,sig,alpha,h);
%               [f,spp,snn]=maternspec(N,sig,alpha,c);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2013--2015 J.M. Lilly --- type 'help jlab_license' for details

if strcmpi(varargin{1}, '--f')
    maternspec_figures,return
end

dt=1;
if length(varargin)==5
    dt=varargin{1};
    varargin=varargin(2:end);
end

N=floor(varargin{1});
omega=fourier(dt,N);
%omega=fourier(N);

sig=varargin{2};
alpha=varargin{3};
H=real(varargin{4});
omegao=-imag(varargin{4});
arrayify(sig,alpha,omegao,H);
%A=frac(A,dt.^(alpha+1/2));
%A=A.*sqrt(dt);

Spp=zeros(length(omega),length(H));
Snn=zeros(length(omega),length(H));

for i=1:length(H)
    if sig(i)~=0
        Spp(:,i)=maternspec_spec(omega,sig(i),omegao(i),H(i),alpha(i));
        if omegao(i)==0
            Snn(:,i)=Spp(:,i);
        else
            Snn(:,i)=maternspec_spec(omega,sig(i),-omegao(i),H(i),alpha(i));
        end
    end
end

if ~isempty(N)
    varargout{1}=omega;
    varargout{2}=Spp;
    varargout{3}=Snn;
else
    varargout{1}=Spp;
    varargout{2}=Snn;
end


function[S]=maternspec_spec(omega,sig,omegao,H,alpha)

d=frac(2*H.^(2*alpha-1),maternc(alpha));
S=frac(sig.^2,((omega-omegao).^2+H.^2).^alpha).*d;


function[]=maternspec_figures
 

%/*************************************************************************
N=1000;

alpha=[1.1 1.5 2 4 8]';
h=1;

figure
subplot(1,3,1)
[f,spp,snn]=maternspec(N,1,alpha,h);
plot([-flipud(f);f]/pi,[flipud(snn);spp])
xlim([-1 1])
title('Matern spectrum')
xlabel('Frequency (cycles/point)')

subplot(1,3,2)
to=[-10:.001:10]';
R=materncov(to,1,alpha,h);
plot(to,R)

title('Matern autocovariance')
xlabel('Time'),xlim([-7 7 ])

subplot(1,3,3)
to=[-2:.001:14]';
g=maternimp(to,alpha,h);
plot(to,g)
title('Matern Green''s function')
xlabel('Time'),xlim([-1.5 12.5 ])


for i=1:3
    subplot(1,3,i)
    ytick off
    %linestyle 2K G k-- K 2G
    vlines(0,'k:')
end
packfig(1,3,'columns')
letterlabels(2)
legend(' \alpha=1.1',' \alpha=1.5',' \alpha=2',' \alpha=4',' \alpha=8')
fontsize 18 14 14 14
set(gcf,'paperposition',[1 1 6 3])
%print -depsc maternspectrum.eps
%\*************************************************************************


%/*************************************************************************
N=1000;

alpha=1;
%alpha=[1.1 1.5 2 4 8]';
%h=1;
h=1./[3 8]';
omegao=1;

figure
subplot(1,3,1)
[f,spp,snn]=maternspec(N,1,alpha,h-sqrt(-1)*omegao);
plot([-flipud(f);f]/pi,[flipud(snn);spp])
xlim([-1 1])
title('Complex OU spectrum')
xlabel('Frequency (cycles/point)')
linestyle k 2G
vlines(1/pi,'D')

subplot(1,3,2)
to=10*[-10:.01:10]';
R=materncov(to,1,alpha,h-sqrt(-1)*omegao);
uvplot(to,R+6*vrep([0 1]*(1+sqrt(-1)),size(R,1),1))
title('Complex OU autocovariance')
xlabel('Time'),xlim([-3.75 3.75 ]*10),ylim([-1.75 10.25])
linestyle k 2G k-- G--

subplot(1,3,3)
to=10*[-2:.01:12]';
g=maternimp(to,alpha,h-sqrt(-1)*omegao);
uvplot(to,4*g+6*vrep([0 1]*(1+sqrt(-1)),size(g,1),1))
title('Complex OU Green''s function')
xlabel('Time'),xlim([-1.75 5.75 ]*10),ylim([-1.75 10.25])
linestyle k 2G k-- G--

for i=1:3
    subplot(1,3,i)
    ytick off
    vlines(0,'k:')
end

packfig(1,3,'columns')
letterlabels(2)
legend(' \omega_o/h=3',' \omega_o/h=8')
fontsize 18 14 14 14
set(gcf,'paperposition',[1 1 6 3])
%print -depsc ouspectrum.eps

%\*************************************************************************
    

%uv=fillbad(uv);
%uv=uv-mean(uv);

% tr=1/10000;
% %tr=1;
% [a,alpha,h] = maternfit(tr,uv);
% [psi,lambda] = sleptap(length(uv),3);
% [f,spp,snn,spn] = mspec(tr,uv,psi,lambda,'adaptive');
% [fm,sppm,snnm] = maternspec(tr,length(uv),a,alpha,h);
% 
% clf
% h=twospecplot(f,[spp sppm],[snn snnm]);
% axes(h(1)),xlin,axes(h(2)),xlin

