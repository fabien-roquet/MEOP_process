function [argo_qc]=ARGO_load_qc(fileIn,varargin);
%function: 	ARGO_load_qc
%object:	read Argo netcdf data files
%
%usage:		[argo_qc]=ARGO_load_qc(fileIn);
%               ---> loads full data set (adjusted+raw data)
%           [argo_qc]=ARGO_load_qc(fileIn,adjusted);
%               ---> adjusted == 0 : load raw data (keep flagged data)
%                                1 : load adjusted data, but keeping headers of all profiles
%                                2 : load raw+adjusted data [default value]
%                                3 : load adjusted data only, removing headers of bad profiles
%                                4 : load raw data only, removing flagged data and headers of bad profiles
%
%           [argo_qc]=ARGO_load_qc(fileIn,adjusted,loadattr);
%               ---> loadattr == 0/[1] : load global attributes
%
%note:      - replaces missing values with NaN
%   		- adds a couple fields: np, nr, list_descr, Pmask, Tmask, Smask
%
%
%inputs:	fileIn		data file name
%
%outputs:	argo_qc	structure containing the various qc fields



% check that file exists and add prefix and suffix if necessary
[pathstr, name, ext] = fileparts(fileIn);
if isempty(pathstr) | strcmp(pathstr,'.'), pathstr=pwd; end
if isempty(ext) | ~strcmp(ext,'.nc'), ext='.nc'; end
fileIn=[pathstr '/' name ext];
if ~exist(fileIn,'file'), error([fileIn ' : file not found']); end

switch nargin,
    case 1,
        adjusted=2;
        loadattr=1;
    case 2,
        adjusted=varargin{1};
        loadattr=1;
    case 3,
        adjusted=varargin{1};
        loadattr=varargin{2};
    otherwise
        error(' ');
end

argo_qc=[];
list_basevar = {'PLATFORM_NUMBER','LATITUDE','LONGITUDE',...
    'JULD_LOCATION','PI_NAME','CYCLE_NUMBER'};

% set list of data to be loaded
switch adjusted
    case 0
        list_var={list_basevar{:},...
            'PRES','TEMP','PSAL',...
            'PRES_QC','TEMP_QC','PSAL_QC'...
            'CHLA','DOXY','LIGHT','CHLA_QC','DOXY_QC','LIGHT_QC'};
    case 1
        list_var={list_basevar{:},...
            'PRES_ADJUSTED','TEMP_ADJUSTED','PSAL_ADJUSTED',...
            'PRES_ADJUSTED_QC','TEMP_ADJUSTED_QC','PSAL_ADJUSTED_QC'...
            'CHLA_ADJUSTED','DOXY_ADJUSTED','LIGHT_ADJUSTED','CHLA_ADJUSTED_QC','DOXY_ADJUSTED_QC','LIGHT_ADJUSTED_QC'};
    case 2,
        list_var={list_basevar{:},...
            'PRES_ADJUSTED','TEMP_ADJUSTED','PSAL_ADJUSTED',...
            'PRES_ADJUSTED_QC','TEMP_ADJUSTED_QC','PSAL_ADJUSTED_QC',...
            'PRES','TEMP','PSAL','TEMP_ADJUSTED_ERROR','PSAL_ADJUSTED_ERROR'...
            'CHLA_ADJUSTED','DOXY_ADJUSTED','LIGHT_ADJUSTED','CHLA_ADJUSTED_QC','DOXY_ADJUSTED_QC','LIGHT_ADJUSTED_QC',...
            'CHLA','DOXY','LIGHT','CHLA_QC','DOXY_QC','LIGHT_QC'};
    case 3
        list_var={list_basevar{:},...
            'PRES_ADJUSTED','TEMP_ADJUSTED','PSAL_ADJUSTED',...
            'PRES_ADJUSTED_QC','TEMP_ADJUSTED_QC','PSAL_ADJUSTED_QC'...
            'CHLA_ADJUSTED','DOXY_ADJUSTED','LIGHT_ADJUSTED','CHLA_ADJUSTED_QC','DOXY_ADJUSTED_QC','LIGHT_ADJUSTED_QC'};
    
    case 4
        list_var={list_basevar{:},...
            'PRES','TEMP','PSAL',...
            'PRES_QC','TEMP_QC','PSAL_QC'...
            'CHLA','DOXY','LIGHT','CHLA_QC','DOXY_QC','LIGHT_QC'};
end
finfo = ncinfo(fileIn);
vars={finfo.Variables.Name};
list_var = intersect(vars,list_var);

% load data
dd=ncload_struct(fileIn,list_var{:});

% handle ADJUSTED fields
if adjusted==1 | adjusted==3
    dd.PRES = dd.PRES_ADJUSTED; dd=rmfield(dd,'PRES_ADJUSTED');
    dd.TEMP = dd.TEMP_ADJUSTED; dd=rmfield(dd,'TEMP_ADJUSTED');
    if isfield( dd,'PSAL_ADJUSTED')
        dd.PSAL = dd.PSAL_ADJUSTED; dd=rmfield(dd,'PSAL_ADJUSTED');
    end
    if isfield( dd,'CHLA_ADJUSTED')
        dd.CHLA = dd.CHLA_ADJUSTED; dd=rmfield(dd,'CHLA_ADJUSTED');
    end
    if isfield (dd,'DOXY_ADJUSTED')
        dd.DOXY = dd.DOXY_ADJUSTED; dd=rmfield(dd,'DOXY_ADJUSTED');
    end
     if isfield (dd,'LIGHT_ADJUSTED')
        dd.LIGHT = dd.LIGHT_ADJUSTED; dd=rmfield(dd,'LIGHT_ADJUSTED');
    end
end

% handle masks
if adjusted==0 || adjusted==4,
    dd.PRES_QC = double(dd.PRES_QC);
    dd.TEMP_QC = double(dd.TEMP_QC);
    if isfield (dd,'PSAL')
        dd.PSAL_QC = double(dd.PSAL_QC);
    end
    if isfield (dd,'CHLA')
        dd.CHLA_QC = double(dd.CHLA_QC);
    end
    if isfield (dd,'DOXY')
        dd.DOXY_QC = double(dd.DOXY_QC);
    end
    if isfield (dd,'LIGHT')
        dd.LIGHT_QC = double(dd.LIGHT_QC);
    end
else
    dd.PRES_QC = double(dd.PRES_ADJUSTED_QC); dd=rmfield(dd,'PRES_ADJUSTED_QC');
    dd.TEMP_QC = double(dd.TEMP_ADJUSTED_QC); dd=rmfield(dd,'TEMP_ADJUSTED_QC');
    if isfield (dd,'PSAL')
        dd.PSAL_QC = double(dd.PSAL_ADJUSTED_QC); dd=rmfield(dd,'PSAL_ADJUSTED_QC');
    end
    if isfield (dd,'CHLA')
        dd.CHLA_QC =double(dd.CHLA_ADJUSTED_QC); dd=rmfield(dd,'CHLA_ADJUSTED_QC');
    end
    if isfield (dd,'DOXY')
        dd.DOXY_QC = double(dd.DOXY_ADJUSTED_QC); dd=rmfield(dd,'DOXY_ADJUSTED_QC');
    end
    if isfield (dd,'LIGHT')
        dd.LIGHT_QC = double(dd.LIGHT_ADJUSTED_QC); dd=rmfield(dd,'LIGHT_ADJUSTED_QC');
    end
end
dd.PRES_QC(dd.PRES_QC>47) = dd.PRES_QC(dd.PRES_QC>47)-48;dd.PRES_QC(dd.PRES_QC==32) = 9;
dd.TEMP_QC(dd.TEMP_QC>47) = dd.TEMP_QC(dd.TEMP_QC>47)-48;dd.TEMP_QC(dd.TEMP_QC==32) = 9;
if isfield (dd,'PSAL')
    dd.PSAL_QC(dd.PSAL_QC>47) = dd.PSAL_QC(dd.PSAL_QC>47)-48;dd.PSAL_QC(dd.PSAL_QC==32) = 9;
end
if isfield (dd,'CHLA')
    dd.CHLA_QC(dd.CHLA_QC>47) = dd.CHLA_QC(dd.CHLA_QC>47)-48;dd.CHLA_QC(dd.CHLA_QC==32) = 9;
end
if isfield (dd,'DOXY')
    dd.DOXY_QC(dd.DOXY_QC>47) = dd.DOXY_QC(dd.DOXY_QC>47)-48;dd.DOXY_QC(dd.PRES_QC==32) = 9;
end
if isfield (dd,'LIGHT')
    dd.LIGHT_QC(dd.LIGHT_QC>47) = dd.LIGHT_QC(dd.LIGHT_QC>47)-48;dd.LIGHT_QC(dd.PRES_QC==32) = 9;
end

dd.Tmask = double(dd.TEMP_QC==1);
if isfield (dd,'PSAL')
    dd.Smask = double(dd.PSAL_QC==1);
else
    dd.Smask=0;
end
% if exist ('dd.CHLA')
%     dd.CHLA = dd.CHLA_ADJUSTED; dd=rmfield(dd,'CHLA_ADJUSTED');
% end
% if exist ('dd.DOXY')
%     dd.DOXY = dd.DOXY_ADJUSTED; dd=rmfield(dd,'DOXY_ADJUSTED');
% end
dd.mask = sum( dd.Tmask+dd.Smask , 1)';
Iprof=1:length(dd.LATITUDE);
if adjusted == 3 || adjusted == 4,
    Iprof=Iprof(find(dd.mask(Iprof)));
end

list_var={'LATITUDE','LONGITUDE','JULD_LOCATION','CYCLE_NUMBER'};
for ii=1:length(list_var),
    if isfield(dd,list_var{ii})
        data=eval(['dd.' list_var{ii} '(Iprof)']);
        argo_qc=setfield(argo_qc,list_var{ii},data);
    end
end

% change date format
argo_qc.JULD_LOCATION=argo_qc.JULD_LOCATION+datenum(1950,1,1);
argo_qc.JULD=argo_qc.JULD_LOCATION;

list_var={'PLATFORM_NUMBER','PI_NAME'};
for ii=1:length(list_var),
    if isfield(dd,list_var{ii})
        data=eval(['dd.' list_var{ii} '(:,Iprof)']);
        argo_qc=setfield(argo_qc,list_var{ii},data);
    end
end

list_var={'PRES','TEMP','PSAL','PRES_QC','TEMP_QC','PSAL_QC','CHLA','DOXY','LIGHT','CHLA_QC','DOXY_QC','LIGHT_QC'};
for ii=1:length(list_var),
    if isfield(dd,list_var{ii})
        data=eval(['dd.' list_var{ii} '(:,Iprof)']);
        argo_qc=setfield(argo_qc,list_var{ii},data);
    end
end

if adjusted==2,
    list_var={'PRES_ADJUSTED','TEMP_ADJUSTED','PSAL_ADJUSTED','TEMP_ADJUSTED_ERROR','PSAL_ADJUSTED_ERROR'...
        ,'CHLA_ADJUSTED','DOXY_ADJUSTED','LIGHT_ADJUSTED','CHLA_ADJUSTED_ERROR','DOXY_ADJUSTED_ERROR','LIGHT_ADJUSTED_ERROR'};
    for ii=1:length(list_var),
        if isfield(dd,list_var{ii})
            data=eval(['dd.' list_var{ii} '(:,Iprof)']);
            argo_qc=setfield(argo_qc,list_var{ii},data);
        end
    end
end

if adjusted~=0,
    argo_qc.PRES(argo_qc.PRES_QC>1)=NaN;
    argo_qc.TEMP(argo_qc.TEMP_QC>1)=NaN;
    if isfield (dd,'PSAL')
        argo_qc.PSAL(argo_qc.PSAL_QC>1)=NaN;
    end
    if isfield (dd,'CHLA')
        argo_qc.CHLA(argo_qc.CHLA_QC>1)=NaN;
    end
    if isfield (dd,'DOXY')
        argo_qc.DOXY(argo_qc.DOXY_QC>1)=NaN;
    end
    if isfield (dd,'LIGHT')
        argo_qc.LIGHT(argo_qc.LIGHT_QC>1)=NaN;
    end

end

if adjusted==2,
    argo_qc.PRES_ADJUSTED(argo_qc.PRES_QC>1)=NaN;
    argo_qc.TEMP_ADJUSTED(argo_qc.TEMP_QC>1)=NaN;
    if exist ('dd.PSAL')
        argo_qc.PSAL_ADJUSTED(argo_qc.PSAL_QC>1)=NaN;
    end
    if exist ('dd.CHLA')
        argo_qc.CHLA_ADJUSTED(argo_qc.CHLA_QC>1)=NaN;
    end
    if exist ('dd.DOXY')
        argo_qc.DOXY_ADJUSTED(argo_qc.DOXY_QC>1)=NaN;
    end
     if exist ('dd.LIGHT')
        argo_qc.LIGHT_ADJUSTED(argo_qc.LIGHT_QC>1)=NaN;
    end
end


%add a couple things:
%--------------------
argo_qc.np=size(argo_qc.PRES_QC,2);
argo_qc.nr=size(argo_qc.PRES_QC,1);
argo_qc.platform_number = cellstr(argo_qc.PLATFORM_NUMBER');
[list_descr,m,n]=unique(argo_qc.platform_number);
argo_qc.list_descr=list_descr;
for kk=1:length(list_descr),
    I=find(strcmp(list_descr{kk},argo_qc.platform_number));
    argo_qc.nprof(kk) = length(I);
    argo_qc.index_tag(I) = kk;
end
argo_qc.ntag=length(list_descr);

if loadattr,
    % load global attributes
    data_att=ncloadatt_struct(fileIn);
    listatt = fieldnames(data_att);
    for ii=1:length(listatt),
        argo_qc = setfield(argo_qc,listatt{ii},getfield(data_att,listatt{ii}));
    end
end



