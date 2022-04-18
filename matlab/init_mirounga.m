function conf = init_mirounga()

conf = init_config;

%% read list deployment
name_list=[conf.processdir 'list_deployment.csv'];
if ~exist(name_list,'file')
    disp('WARNING: the file list_deployment.csv was not found!')
    disp(['path: ' name_list])
    return
end
list_deployment = readtable(name_list,'ReadRowNames',1,'Delimiter',',');
list_deployment_code = list_deployment.deployment_code;

% remove old repeated deployments
list_old_repeated_deployments = {'ct10r','ct12r','ct17r'};
for depl = list_old_repeated_deployments,
    if any(strcmp(depl{1},list_deployment_code)),
        list_deployment(depl,:)=[];
    end
end

conf.list_deployment = list_deployment;
conf.list_deployment_code = list_deployment.deployment_code;

%% pre-process old deployments (done once and for all)
if 1==0,
    preprocessing_old_deployments;
end

%% load deployment information
name_json = [conf.json 'deployment3.json'];
filetext = fileread(name_json);
deployment_json = jsondecode(filetext);

name_json = [conf.json 'deployment2_patch.json'];
filetext = fileread(name_json);
deployment_patch_json = jsondecode(filetext);
deployment_code_json={deployment_json(:).deployment_code};
deployment_code_patch_json={deployment_patch_json(:).deployment_code};
for kk=1:length(deployment_code_patch_json),
    if ~ismember(deployment_code_patch_json{kk},deployment_code_json)
        deployment_json=[deployment_json;deployment_patch_json(kk)];
    end
end

% add new deployments
new_deployments = 0;
for ii=1:length(deployment_json)
    if ~any(strcmp(deployment_json(ii).deployment_code,conf.list_deployment.deployment_code)),
        % list_deployment.csv :
        % deployment_code,pi_code,process,public,country,task_done,
        % first_version,last_version,start_date,end_date,start_date_jul
        EXP = deployment_json(ii).deployment_code;
        if any(strcmp(EXP,list_old_repeated_deployments)),
            continue
        end
        conf.list_deployment(EXP,:)={'',1,0,'UNKNOWN','','',NaT,NaT,NaN};
        conf.list_deployment{EXP,'pi_code'} = {deployment_json(ii).pi_code};
        disp(sprintf('%s added in list_deployment.csv',deployment_json(ii).deployment_code));
        new_deployments = 1;
    end
end

conf.deployment_json = deployment_json;
conf.list_deployment_json = {deployment_json(:).deployment_code};


%% load platform information
name_json = [conf.json 'platform3.json'];
filetext = fileread(name_json);
platform_json = jsondecode(filetext);
deployment_code_json=unique({platform_json(:).deployment_code});

name_json = [conf.json 'platform2_patch.json'];
filetext = fileread(name_json);
platform_patch_json = jsondecode(filetext);
for kk=1:length(platform_patch_json),
    platform_patch_json(kk).cond_sensor='';
    if ismember(platform_patch_json(kk).deployment_code,conf.list_deployment_code) ...
            && ~ismember(platform_patch_json(kk).deployment_code,deployment_code_json),
        platform_json=[platform_json;platform_patch_json(kk)];
    end
end

conf.platform_json = platform_json;
conf.list_smru_platform_code = {conf.platform_json(:).smru_platform_code};


%% add platform information based on deployment.json
for kk=1:length(conf.platform_json),
    deployment_code = conf.platform_json(kk).deployment_code;
    I = find(strcmp(deployment_code,conf.list_deployment_json));
    if I,
        conf.platform_json(kk).description = conf.deployment_json(I).description;
        conf.platform_json(kk).gts = conf.deployment_json(I).gts;
    else
        conf.platform_json(kk).description ='';
        conf.platform_json(kk).gts ='';
    end
end


%% build group list
conf.list_group=unique(conf.list_deployment.country);
conf.list_color= linspecer(length(conf.list_group));
conf.Itag_group={};
for kk=1:length(conf.list_group),
    conf.Itag_group{kk}=find(strcmp(conf.list_group{kk},conf.list_deployment.country));
end


%% conf.start_date for each deployment
aux = {conf.platform_json.deployment_code;
    conf.platform_json.time_coverage_start;
    conf.platform_json.time_coverage_end};
for kEXP = 1:length(conf.list_deployment_code),
    EXP = conf.list_deployment_code{kEXP};
    I = find(strcmp(aux(1,:),EXP));
    if I & isempty(conf.list_deployment{EXP,'end_date'}),
        if aux{2,I(1)}, conf.list_deployment{EXP,'start_date'} = {aux{2,I(1)}(1:10)}; end
        if aux{3,I(1)}, conf.list_deployment{EXP,'end_date'}   = {aux{3,I(1)}(1:10)}; end
    end
end

% trim pi_code
pi_code = conf.list_deployment.pi_code;
for kk = 1:length(pi_code),
    pi_code{kk} = strtrim(pi_code{kk});
end
conf.list_deployment.pi_code = pi_code;

% add jul_date of start_date
date = conf.list_deployment.start_date;
juldate = zeros(size(date));
for kk = 1:length(juldate),
    if ~isempty(date(kk)), juldate(kk) = datenum(date(kk)); end
end
conf.list_deployment.start_date_jul = juldate;

% update list_deployment.json
name_list=[conf.processdir 'list_deployment.csv'];
writetable(conf.list_deployment,name_list,'WriteRowNames',1,'Delimiter',',');
if new_deployments | ismember(conf.list_group,'UNKNOWN'),
    error(['Update information for new deployments in ' name_list]);
end

%% read list deployment CTD HR
name_list=[conf.processdir 'list_deployment_hr.csv'];
if ~exist(name_list,'file')
    error(['WARNING: the file ' name_list ' was not found!']);
end
conf.list_deployment_hr = readtable(name_list,'ReadRowNames',1,'Delimiter',',');
conf.hr_smru_name = conf.list_deployment_hr.smru_platform_code;
conf.hr_continuous = conf.list_deployment_hr.continuous;

%% load config files
name_file=[conf.processdir 'table_coeff.csv'];
if ~exist(name_file,'file')
    error(['WARNING: the file ' name_file ' was not found!'])
end
conf.table_coeff = readtable(name_file,'ReadRowNames',1,'Delimiter',',');

name_file=[conf.processdir 'table_param.csv'];
if ~exist(name_file,'file')
    error(['WARNING: the file ' name_file ' was not found!'])
end
conf.table_param = readtable(name_file,'ReadRowNames',1,'Delimiter',',');

name_file=[conf.processdir 'table_salinity_offsets.csv'];
if ~exist(name_file,'file')
    error(['WARNING: the file ' name_file ' was not found!'])
end
conf.table_salinity_offsets = readtable(name_file,'ReadRowNames',1,'Delimiter',',');

name_file=[conf.processdir 'table_meta.csv'];
if ~exist(name_file,'file')
    error(['WARNING: the file ' name_file ' was not found!'])
end
conf.table_meta = readtable(name_file,'ReadRowNames',1,'Delimiter',',');

name_file=[conf.processdir 'table_filter.csv'];
if ~exist(name_file,'file')
    error(['WARNING: the file ' name_file ' was not found!'])
end
conf.table_filter = readtable(name_file,'ReadRowNames',0,'Delimiter',',');

name_file=[conf.processdir 'table_split_tags.csv'];
if ~exist(name_file,'file')
    error(['WARNING: the file ' name_file ' was not found!'])
end
conf.table_split_tags = readtable(name_file,'ReadRowNames',1,'Delimiter',',');


%% load information about updated loc
conf.crawl.list    = dir([conf.crawl.locdir,'*_argos_crawl.csv']);
conf.crawl.ptt     = {};
conf.crawl.smru_name     = {};
for kk=1:length(conf.crawl.list)
    c = strsplit(conf.crawl.list(kk).name,{'_','.'});
    conf.crawl.smru_name{kk} = c{2};
    conf.crawl.ptt{kk} = c{3};
end

conf.cls.list    = dir([conf.cls.locdir,'*.smoothing.csv']);
conf.cls.ptt     = {};
conf.cls.datemin = {};
conf.cls.datemax = {};
for kk=1:length(conf.cls.list)
    c = strsplit(conf.cls.list(kk).name,{'_','.'});
    conf.cls.ptt{kk} = c{1};
    conf.cls.datemin{kk} = c{2};
    conf.cls.datemax{kk} = c{3};
end
conf.cls.datemin_jul = datenum(conf.cls.datemin);
conf.cls.datemax_jul = datenum(conf.cls.datemax);


end




