function[meandt,sigdt,meddt,maxdt]=sampletimes(num)
%SAMPLETIMES  Computes mean sampling intervals and their statistics.
%
%   DT=SAMPLETIMES(NUM) where NUM is a column vector of times, returns the 
%   mean interval between sample times DT.
%
%   [DT,SIG,MED,MAX]=SAMPLETIMES(NUM) also returns the standard deviation
%   of sample intervals SIG, the median sample interval MED, and the
%   maximum sample interval MAX.
%
%   IF NUM is a cell array of numerical arrays, the output of SAMPLETIMES
%   will be arrays of size LENGTH(NUM) x 1.
%
%   Usage: dt=sampletimes(num); 
%          [dt,sig,med,max]=sampletimes(num);
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014 J.M. Lilly --- type 'help jlab_license' for details
 
% if strcmpi(num, '--t')
%     sampletimes_test,return
% end

if ~iscell(num)
    [meandt,sigdt,meddt]=sampletimes_one(num);
else
    [meandt,sigdt,meddt]=vzeros(length(num),1);
    for i=1:length(num)
        if nargout>3
            [meandt(i),sigdt(i),meddt(i),maxdt(i)]=sampletimes_one(num{i});
        elseif nargout>2
            [meandt(i),sigdt(i),meddt(i)]=sampletimes_one(num{i});
         elseif nargout>1
            [meandt(i),sigdt(i)]=sampletimes_one(num{i});
         else
            meandt(i)=sampletimes_one(num{i});
         end
    end
end

function[meandt,sigdt,meddt,maxdt]=sampletimes_one(num)

if isempty(num)
   meandt=nan;
   sigdt=nan;
   meddt=nan;
else
    dt=num(2:end,:)-num(1:end-1,:);
    
    meandt=vmean(dt,1);
    if nargout>1
        sigdt=vstd(dt,1);
    end
    if nargout >2
        meddt=vmedian(dt,1);
    end
    if nargout >3
        maxdt=maxmax(dt);
    end
end
%function[]=sampletime_test
 
%reporttest('SAMPLETIME',aresame())
