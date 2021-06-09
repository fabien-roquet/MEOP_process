function info_deployment=load_info_deployment(conf,EXP,one_smru_name)

%% load metadata deployment :
% EXP, PI, ODV name and directory name 

if isempty(conf),
    conf = init_mirounga;
end

if ~exist('one_smru_name','var') % all tags from EXP deployment
    one_smru_name = '';
elseif isempty(EXP),
    EXP=EXP_from_smru_name(one_smru_name);
end

if iscell(EXP) && length(EXP)==1,
    EXP = EXP{1};
end

if ~ismember(EXP,conf.list_deployment_code)
% deployment_code,pi_code,process,public,country,task_done,first_version,last_version,start_date,end_date,start_date_jul
% {"deployment_code":"hg63","description":"Hammill Hg CTD 2019","pi_code":"HAMMILL ","gts":"Y","dt_created":"2019-02-14T10:07:07Z","dt_modified":"2020-01-08T16:03:27Z"}
    info_deployment.EXP      = EXP;
    info_deployment.list_tag = [];
    info_deployment.invalid_code = 1;
    disp([EXP ' is not a valid deployment code. Update list_deployment.csv']);
    return
end

info_deployment.EXP      = EXP;
info_deployment.PI       = conf.list_deployment{EXP,'pi_code'}{1};
info_deployment.PI       = info_deployment.PI(info_deployment.PI~=' ');
info_deployment.NATION   = conf.list_deployment{EXP,'country'}{1};
info_deployment.nomfic   = [EXP '_ODV.txt'];
if ismember(EXP,{'ct3','ct7','ct11','wd3'}),  
    info_deployment.nomfic   = [EXP '_fcell.mat']; 
end
info_deployment.dir      = [conf.datadir,EXP,'/'];

info_deployment.process   = conf.list_deployment{EXP,'process'};
info_deployment.public    = conf.list_deployment{EXP,'public'};


info_deployment.list_tag = dir([info_deployment.dir one_smru_name '*_lr0_prof.nc']);
info_deployment.list_tag_lr0 =dir([info_deployment.dir one_smru_name '*_lr0_prof.nc']);
info_deployment.list_tag_lr1 =dir([info_deployment.dir one_smru_name '*_lr1_prof.nc']);
info_deployment.list_tag_hr0 =dir([info_deployment.dir one_smru_name '*_hr0_prof.nc']);
info_deployment.list_tag_hr1 =dir([info_deployment.dir one_smru_name '*_hr1_prof.nc']);
info_deployment.list_tag_hr2 =dir([info_deployment.dir one_smru_name '*_hr2_prof.nc']);
info_deployment.list_tag_fr0 =dir([info_deployment.dir one_smru_name '*_fr0_prof.nc']);
info_deployment.list_tag_fr1 =dir([info_deployment.dir one_smru_name '*_fr1_prof.nc']);

info_deployment.list_smru_name = {info_deployment.list_tag.name};
count = 0;
for a = info_deployment.list_smru_name,
    count = count+1;
    info_deployment.list_smru_name{count} = a{1}(1:end-12);
end


