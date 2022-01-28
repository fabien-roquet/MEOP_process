function create_ncargo(conf,EXP,one_smru_name)


if isempty(conf),
    conf = init_mirounga;
end

if ~exist('one_smru_name','var')  % all tags from EXP deployment
    one_smru_name = '';
elseif isempty(EXP),
    EXP=EXP_from_smru_name(one_smru_name);
end


% don't process it if no raw odv file
info_deployment=load_info_deployment(conf,EXP);
if ~exist([conf.rawdir info_deployment.nomfic]),
    error(sprintf('%s: no raw file. not processed.',EXP));
end

%  create directory
[s,mess,messid] = mkdir(info_deployment.dir);

%  open diary file
diary_file = [info_deployment.dir EXP '_diary.txt'];
if ~exist('one_smru_name','var') | isempty(one_smru_name) % all tags from EXP deployment
    one_smru_name = '';
    if exist(diary_file,'file'), delete(diary_file); end
    diary(diary_file)
    disp(['Process EXP=' EXP])
else  % tag smru_tag only
    info_deployment=load_info_deployment(conf,EXP,one_smru_name);
    diary(diary_file)
    disp(['Process smru_name=' one_smru_name])
end


%  create ncARGO file lr0
if ismember(info_deployment.EXP,{'ct3','ct7','ct11','wd3'})
    sc_load_KERold2prof;
else
    sc_load_odv2prof;
end


%  write global attributes
info_deployment=load_info_deployment(conf,EXP,one_smru_name);
for ktag=1:length(info_deployment.list_smru_name),
    smru_name = info_deployment.list_smru_name{ktag};
    name_prof = sprintf('%s%s_lr0_prof.nc',info_deployment.dir,smru_name);
    sc_write_global_attribute;
end


% create default coefficients if needed
info_deployment=load_info_deployment(conf,EXP,one_smru_name);
list_tag = info_deployment.list_smru_name;
for ktag=1:length(list_tag),
    if ~ismember(conf.table_coeff.Properties.RowNames,list_tag{ktag}),
        new_tag = {0,0,0,0,0,0,'no comment'};
        conf.table_coeff=[conf.table_coeff;new_tag];
        conf.table_coeff.Properties.RowNames{end}=list_tag{ktag};
        name_file=[conf.processdir 'table_coeff.csv'];
        writetable(conf.table_coeff,name_file,'WriteRowNames',1,'Delimiter',',');
    end
end


%  apply filters
info_deployment=load_info_deployment(conf,EXP,one_smru_name);
for ktag=1:length(info_deployment.list_smru_name),    
    smru_name = info_deployment.list_smru_name{ktag};
    suffix = '_lr0';
    name_prof = sprintf('%s%s%s_prof.nc',info_deployment.dir,smru_name,suffix);
    sc_filtre_seals_qc;
end


%  adjust locations
sc_adjust_locations;


%  create ncARGO file hr0
sc_create_hr0;


diary off

