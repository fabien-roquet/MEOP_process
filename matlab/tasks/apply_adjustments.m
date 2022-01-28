function apply_adjustments(conf,EXP,one_smru_name)


if isempty(conf),
    conf = init_mirounga;
end

if ~exist('one_smru_name','var') % all tags from EXP deployment
    one_smru_name = '';
elseif isempty(EXP),
    EXP=EXP_from_smru_name(one_smru_name);
end


info_deployment=load_info_deployment(conf,EXP,one_smru_name);
if ~exist(info_deployment.dir), return, end

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

% apply offset
for kk=1:length(info_deployment.list_tag),
    
    name_lr0 = sprintf('%s%s',info_deployment.dir,info_deployment.list_tag(kk).name);
    name_lr1 = strrep(name_lr0, '_lr0', '_lr1' );
    name_hr0 = strrep(name_lr0, '_lr0', '_hr0' );
    name_hr1 = strrep(name_lr0, '_lr0', '_hr1' );
    name_fr0 = strrep(name_lr0, '_lr0', '_fr0' );
    name_fr1 = strrep(name_lr0, '_lr0', '_fr1' );
    aux = {name_lr0,name_lr1,name_hr0,name_hr1,name_fr0,name_fr1};
    
    list_files={};
    for ii = 1:length(aux),
        if exist(aux{ii},'file'),
            list_files = { list_files{:}, aux{ii}};
        end
    end
    
    if isempty(list_files),
        continue;
    end
    
    smru_name = info_deployment.list_smru_name{kk};
    for ii = 1:length(list_files),
        
        name_prof=list_files{ii};
        
        coeff = conf.table_coeff;
        salinity_offsets = conf.table_salinity_offsets;
        correct_offset_qc(name_prof,smru_name,coeff,salinity_offsets);
        
        temp_error = conf.table_param{EXP,'temp_error'};
        psal_error = conf.table_param{EXP,'psal_error'};
        assign_error_estimate(name_prof,name_lr0,temp_error,psal_error)
        
    end
    
end


