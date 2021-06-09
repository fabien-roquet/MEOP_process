function[bool,rq,w,om]=isridgepoint(w,fs,chi,str,fmin,fmax,mask)
%ISRIDGEPOINT  Finds wavelet ridge points using one of several criterion.
%
%   ISRIDGEPOINT is a low-level function called by RIDGEWALK.
%  
%   BOOL=ISRIDGEPOINT(W,FS,CHI,STR) where W is a wavelet transform matrix  
%   at *radian* frequecies FS, finds all ridge points of W with amplitudes
%   |W| exceeding the amplitude cutoff A.  Several different different 
%   ridge defintions may be used and are specified by STR.
%
%   BOOL is a matrix of the same size as W, which is equal to one for 
%   those elements of W which are ridge points, and zero otherwise.
%
%   STR may be either of the following:
%
%        'phase'       Rate of transform change of phase definition
%        'amplitude'   Maxima of transfom amplitude definition
%
%   For all definitions, ISRIDGEPOINT rejects spurious ridge points.
%   These tend to occur on the flanks of interesting signals, and 
%   reflect the wavelet structure rather than the signal structure.
%
%   A ridge point is considered spurious if either it is located at an
%   amplitude minima, or if the frequency anomaly (transform frequency
%   minus scale frequency) is a maximum.
%
%   See also RIDGEQUANTITY, RIDGEWALK.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2015 J.M. Lilly --- type 'help jlab_license' for details
 
%        'groove'      Joint amplitude / phase definition

disp('RIDGEWALK looking for ridge points...')

%Re-doing way of handling amplitude ridges to remove ridge breaking
if size(w,3)~=1
    [a,om]=instmom(w,1,3);
else
    [a,om]=instmom(w);
end

if strcmpi(str(1:3),'amp')
    rq=a;
elseif strcmpi(str(1:3),'pha')
    fsmat=vrep(fs(:)',size(w,1),1);
    rq=om-fsmat;
end    
if size(w,3)~=1
    phaseavg=frac(sum(abs(w).*w,3),sum(abs(w).^2,3));
    w=sqrt(sum(abs(w).^2,3)).*rot(angle(phaseavg));
end

%figure,jpcolor(rq'),shading flat

rqm=circshift(rq,+1,2);
rqp=circshift(rq,-1,2);

if strcmpi(str(1:3),'amp')
   bool=(rqm<=rq)&(rqp<=rq);
elseif strcmpi(str(1:3),'pha')
   %This is d/ds < 0 since scale decreases in columns
   bool=(rqm<0&rqp>=0)|(rqm<=0&rqp>0);
end
%[ii,jj]=find(bool);
%figure,plot(ii,jj,'.')

err=abs(rq);

%Ensure maximum not minimum
bool((bool&circshift(bool,-1,2))&err>circshift(err,-1,2))=0; 
bool((bool&circshift(bool,+1,2))&err>circshift(err,+1,2))=0; 

bool1= ~isnan(w);   %Remove NANs
bool2=~(abs(w)<chi);  %Remove those less than cutoff amplitude
bool=bool.*bool1.*bool2;
bool(:,[1 end])=0;

%Running FMIN and FMAX 
if ~isempty(fmin)
    %figure,plot(fsmat)
    %figure,plot(om)
    %    figure,plot(om)fmin,fmax
    %    hlines(fmin),hlines(fmax)
        
    fmin=vrep(fmin,size(w,2),2);
    fmax=vrep(fmax,size(w,2),2);
    if size(fmin,1)==1
        fmin=vrep(fmin,size(w,1),1);
        fmax=vrep(fmax,size(w,1),1);
    end
    %vsize(fmin,fmax,fsmat)
    
    %This way is more intuitive, and seems better at rejecting short
    %ridges.  Otherwise the results can be identical.
    bool3=(om>fmin)&(om<fmax);
    %bool3=(fsmat>=fmin)&(fsmat<=fmax)&(om>fmin)&(om<fmax);
    %fsmat=vrep(fs(:)',size(w,1),1);
    %bool3=(fsmat>=fmin)&(fsmat<=fmax);
    bool=bool.*bool3;
    %figure,plot(om),hlines(fs)
end

if ~isempty(mask)
    bool=bool.*mask;
end
disp(['RIDGEWALK found ' int2str(length(find(bool))) ' ridge points.'])
