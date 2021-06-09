function [argo_data]=ARGO_load(fileIn,varargin);
%function: 	ARGO_load
%object:	read Argo netcdf data files in the "Argo netcdf format"
%
%usage:		[argo_data]=ARGO_load(fileIn);
%               ---> loads full data set
%           [argo_data]=ARGO_load(fileIn,list_vars);
%               ---> loads additional variables listed in list_vars cell
%               array (e.g. list_vars={'prof_T','prof_Tweight'})
%           [argo_data]=ARGO_load(fileIn,list_vars,list_platform_numbers);
%               ---> loads only data having a PLATFORM_NUMBERs listed in 
%               list_platform_numbers cell array 
%
%note:      - replaces missing values with NaN
%   		- adds a couple fields: np, nr, list_descr, Pmask, Tmask, Smask
%       
%
%inputs:	fileIn		data file name
%
%outputs:	argo_data	structure containing the various fields/vectors



% check that file exists and add prefix and suffix if necessary

[pathstr, name, ext] = fileparts(fileIn);
if isempty(pathstr) | strcmp(pathstr,'.'), pathstr=pwd; end
if isempty(ext) | ~strcmp(ext,'.nc'), ext='.nc'; end
fileIn=[pathstr '/' name ext];
if ~exist(fileIn,'file'), error([fileIn ' : file not found']); end

list_additional_vars={};
%if nargin>1; list_additional_vars=varargin{1}; end

list_var={'PLATFORM_NUMBER','JULD','LATITUDE','LONGITUDE','DIRECTION',...
    'PRES','TEMP','PSAL'};
dd=ncload_struct(fileIn,list_var{:});

Iprof=1:length(dd.JULD);
if nargin>2; 
    list_platform_numbers=varargin{2};
    Iprof=find(ismember(cellstr(dd.PLATFORM_NUMBER'),list_platform_numbers));
end

argo_data=[];

list_var={'PLATFORM_NUMBER'};
for ii=1:length(list_var),
    data=eval(['dd.' list_var{ii} '(:,Iprof)']);
    argo_data=setfield(argo_data,list_var{ii},data);
end

list_var={'JULD','LATITUDE','LONGITUDE','DIRECTION'};
for ii=1:length(list_var),
    data=eval(['dd.' list_var{ii} '(Iprof)']);
    argo_data=setfield(argo_data,list_var{ii},data);
end

list_var={'PRES','TEMP','PSAL'};
for ii=1:length(list_var),
    data=double(eval(['dd.' list_var{ii} '(:,Iprof)']));
    argo_data=setfield(argo_data,list_var{ii},data);
end

if ~isempty(list_additional_vars),
    list_var=list_additional_vars;
    for ii=1:length(list_var),
        ncload(fileIn,list_var{ii});
        if isempty(eval(list_var{ii})), 
            %disp(['no variable named: ' list_var{ii}]);
            continue
        end
        data=eval(list_var{ii});
        if size(data,1)==1,
            data=data(Iprof);
        else
            data=data(:,Iprof);
        end
        argo_data=setfield(argo_data,list_var{ii},data);
    end
end

% change date format
argo_data.JULD=argo_data.JULD+datenum(1950,1,1);

%add a couple things:
%--------------------
argo_data.np=size(argo_data.PRES,2);
argo_data.nr=size(argo_data.PRES,1);
argo_data.platform_number = cellstr(argo_data.PLATFORM_NUMBER');
[list_descr,m,n] = unique(cellstr(argo_data.platform_number'),'stable');
argo_data.ntag      =length(list_descr);
argo_data.list_descr=list_descr;
argo_data.index_tag = n;
for kk=1:length(list_descr),
    argo_data.nprof(kk) = length(find(n==kk));
end


if isfield(argo_data,'SMRU_NAME')
    
    argo_data.smru_name = cellstr(argo_data.SMRU_NAME');
    argo_data.list_descr_smru=argo_data.list_descr;
    for kk=1:length(argo_data.list_descr),
        I=find(strcmp(argo_data.list_descr{kk},argo_data.platform_number));
        if isempty(I), continue; end
        argo_data.list_descr_smru{kk}=argo_data.smru_name{I(1)};
    end
    
end

