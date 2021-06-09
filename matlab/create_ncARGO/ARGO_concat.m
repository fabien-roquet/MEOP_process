function [argo_data]=ARGO_concat(argo_data1,argo_data2);
%   [argo_data]=ARGO_concat(argo_data1,argo_data2);;
%   concatenates argo_data1 and argo_data2
%   Only fields present in both datasets are concatenated.
%
% if argo_data1 is empty, argo_data2 is simply returned

if isempty(argo_data1) || length(argo_data1.LATITUDE)==0,
    argo_data=argo_data2;
elseif isempty(argo_data2) || length(argo_data2.LATITUDE)==0,
    argo_data=argo_data1;
else
    
    fldNames1=fieldnames(argo_data1);np1=length(argo_data1.JULD);
    fldNames2=fieldnames(argo_data2);np2=length(argo_data2.JULD);
    
    fldNames1=setdiff(fldNames1,{'np','nr','nprof','index_tag','ntag','list_descr','Pmask','Tmask','Smask'},'stable');
    fldNames2=setdiff(fldNames2,{'np','nr','nprof','index_tag','ntag','list_descr','Pmask','Tmask','Smask'},'stable');
    
    argo_data=[];
    
    %------------
    for iFld=1:length(fldNames1);
        eval(['tmp1=argo_data1.' fldNames1{iFld} ';']);
        if isfield(argo_data2,fldNames1{iFld})
            eval(['tmp2=argo_data2.' fldNames1{iFld} ';']);
            if isempty(tmp1), tmp1=' '; end
            if isempty(tmp1), tmp2=' '; end
            if iscellstr(tmp1(1))
                argo_data=setfield(argo_data,fldNames1{iFld},{tmp1{:} tmp2{:}}');
            elseif isnumeric(tmp1(1))
                if isvector(tmp1) & argo_data1.np>1 & argo_data2.np>1
                    tmp=zeros(np1+np2,1)*NaN;
                    tmp(1:np1)=tmp1;
                    if length(tmp2)==np2, tmp(np1+1:np1+np2)=tmp2;
                    else, tmp(np1+1:np1+np2)=NaN; end
                    argo_data=setfield(argo_data,fldNames1{iFld},tmp);
                else
                    try
                        size1=max(size(tmp1,1),size(tmp2,1));
                        tmp=zeros(size1,np1+np2)*NaN;
                        tmp(1:size(tmp1,1),1:np1)=tmp1;
                        tmp(1:size(tmp2,1),np1+1:np1+np2)=tmp2;
                        argo_data=setfield(argo_data,fldNames1{iFld},tmp);
                    catch
                        tmp=zeros(np1+np2,1)*NaN;
                        try; tmp(1:np1)=tmp1; catch
                            1
                        end
                        if isempty(tmp2), tmp2=NaN*zeros(np2,1); end
                        tmp(np1+1:np1+np2)=tmp2;
                        argo_data=setfield(argo_data,fldNames1{iFld},tmp);
                    end
                end
            elseif ischar(tmp1(1))
                if strcmp(fldNames1{iFld},'PLATFORM_NUMBER') | strcmp(fldNames1{iFld},'PI_NAME')
                    tmp = [tmp1 tmp2];
                    argo_data=setfield(argo_data,fldNames1{iFld},tmp);
                end
            end
        end
    end
    
end

%add a couple things:
%--------------------
argo_data.np=length(argo_data.JULD);
argo_data.nr=size(argo_data.PRES,1);
[list_descr,m,n]=unique(cellstr(argo_data.platform_number'),'stable');
argo_data.ntag      =length(list_descr);
argo_data.list_descr=list_descr;
argo_data.index_tag = n;
for kk=1:length(list_descr),
    argo_data.nprof(kk) = length(find(n==kk));
end

if isfield(argo_data,'PRES_ADJUSTED')
    argo_data.Pmask=double(~isnan(argo_data.PRES_ADJUSTED));
    argo_data.Tmask=double(~isnan(argo_data.TEMP_ADJUSTED));
    argo_data.Smask=double(~isnan(argo_data.PSAL_ADJUSTED));
else
    argo_data.Pmask=double(~isnan(argo_data.PRES));
    argo_data.Tmask=double(~isnan(argo_data.TEMP));
    argo_data.Smask=double(~isnan(argo_data.PSAL));
end
