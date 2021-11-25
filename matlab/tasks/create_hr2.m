function create_hr2(conf,EXP,one_smru_name)

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

disp(['Create hr2 version [by default, hr1 --> hr2]'])
list_tag = info_deployment.list_tag_hr1;
list_tag_fr1 = info_deployment.list_tag_fr1;
for ll=1:length(list_tag),
    
    name_prof = list_tag(ll).name;
    [smru_prefix,Nsplit,suffix,namedir] = smru_name_from_name_prof(name_prof);
    smru_name = gen_smru_name(smru_prefix,Nsplit);
    name_prof_hr1 = gen_name_prof(smru_prefix,Nsplit,'hr1',info_deployment.dir);
    name_prof_fr1 = gen_name_prof(smru_prefix,Nsplit,'fr1',info_deployment.dir);
    name_prof_hr2 = gen_name_prof(smru_prefix,Nsplit,'hr2',info_deployment.dir);
    
    if ~ismember(smru_name,conf.hr_smru_name) || ~isfile(name_prof_fr1),
        disp(['  ' smru_name ': hr1 --> hr2'])
        copyfile(name_prof_hr1, name_prof_hr2,'f');
        continue;
    end
    
    
    kk = find(strcmp(smru_name,conf.hr_smru_name));
    
    if conf.hr_continuous(kk),
        
        copyfile(name_prof_hr1, name_prof_hr2,'f');
        data_hr2 = ncload_struct(name_prof_hr2);
        data_fr1 = ncload_struct(name_prof_fr1);
        I = round(interp1(data_fr1.JULD,1:length(data_fr1.JULD),data_hr2.JULD));
        J = find(~isnan(I));
        
        list_var1={'JULD','JULD_LOCATION','LATITUDE','LATITUDE'};
        for jj=1:length(list_var1),
            aux  = getfield(data_hr2,list_var1{jj});
            aux1 = getfield(data_fr1,list_var1{jj});
            aux(J) = aux1(I(J));
            ncwrite(name_prof_hr2,list_var1{jj},aux);
        end
        
        list_var2={'PRES','PRES_QC','PRES_ADJUSTED','PRES_ADJUSTED_ERROR',...
            'TEMP','TEMP_QC','TEMP_ADJUSTED','TEMP_ADJUSTED_ERROR',...
            'PSAL','PSAL_QC','PSAL_ADJUSTED','PSAL_ADJUSTED_ERROR'};
        for jj=1:length(list_var2),
            aux  = getfield(data_hr2,list_var2{jj});
            aux1 = getfield(data_fr1,list_var2{jj});
            aux(:,J) = aux1(:,I(J));
            ncwrite(name_prof_hr2,list_var2{jj},aux);
        end
        
        list_var3={'PARAMETER','SCIENTIFIC_CALIB_EQUATION','SCIENTIFIC_CALIB_COEFFICIENT'};
        for jj=1:length(list_var3),
            aux  = getfield(data_hr2,list_var3{jj});
            aux1 = getfield(data_fr1,list_var3{jj});
            aux(:,1:3,:,J) = aux1(:,1:3,:,I(J));
            ncwrite(name_prof_hr2,list_var3{jj},aux);
        end

        disp(['  ' smru_name ': continuous mode fr1 taken at hr1 datetime  --> hr2'])
        
    else
        
        disp(['  ' smru_name ': high resolution fr1 --> hr2'])
        copyfile(name_prof_fr1, name_prof_hr2,'f');
        
    end
    
    
end

