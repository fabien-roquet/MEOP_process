%% copy final data files in public folder

folder_output_SMS = sprintf('%s%s/',conf.public,conf.version_SMS);
[status,msg,msgID] = rmdir(folder_output_SMS, 's');
[s,mess,messid] = mkdir(folder_output_SMS);


[status,message] = copyfile('README_licenseODbl.txt',folder_output_SMS,'f');
[status,message] = copyfile('info_tags_sms.csv',folder_output_SMS,'f');
[status,message] = copyfile(sprintf('%sgroups/FRANCE_SMS_mapSH.png',conf.mapsdir),[folder_output_SMS 'map_MEOP-SMS.png'],'f');
[status,message] = copyfile(sprintf('%sglobal/map_SMS.png',conf.mapsdir),folder_output_SMS,'f');

EXP_all=tags_processed(conf);
list_NATION=unique(EXP_all.country);

for  kNATION=1:length(list_NATION),
    
    NATION = list_NATION{kNATION};
    EXPs=tags_processed(conf,NATION);
    list_EXP = EXPs.deployment_code;
    
    if ~any(EXPs.process==1), continue; end
    if ~any(EXPs.public ==1), continue; end
    
    for kEXP = 1:length(list_EXP),
        
        EXP = list_EXP{kEXP};
        info_deployment=load_info_deployment(conf,EXP);
        list_tag_fr1 = info_deployment.list_tag_fr1;
        Ntag = length(list_tag_fr1);
        for jj=1:Ntag,
            
            ncfile_fr1 = sprintf('%s%s',info_deployment.dir,list_tag_fr1(jj).name);
            smru_name = list_tag_fr1(jj).name(1:end-12);
            exist_continuous = conf.hr_continuous(find(strcmp(smru_name,conf.hr_smru_name)));
            trajfile = sprintf('%s%s_traj.nc',info_deployment.dir,smru_name);
            
            if ~exist_continuous | ~exist(ncfile_fr1,'file'), continue; end
            
            [s,mess,messid] = mkdir([folder_output_SMS NATION '/ncARGO_traj/']);
            [s,mess,messid] = mkdir([folder_output_SMS NATION '/ncARGO_prof/']);
            [s,mess,messid] = mkdir([folder_output_SMS NATION '/MAPS/']);
            [s,mess,messid] = mkdir([folder_output_SMS NATION '/PDF/']);
            
            [status,message] = copyfile('README_licenseODbl.txt',[folder_output_SMS NATION],'f');
            [status,message] = copyfile( ...
                sprintf('%s%s_fr1_doc_adj.pdf',conf.texdir,EXP),...
                [folder_output_SMS NATION '/PDF/'],'f');
            [status,message] = copyfile(ncfile_fr1,[folder_output_SMS NATION '/ncARGO_prof/' smru_name '_SMS_prof.nc'],'f');
            [status,message] = copyfile(trajfile,[folder_output_SMS NATION '/ncARGO_traj/' smru_name '_SMS_traj.nc'],'f');
            [status,message] = copyfile( ...
                sprintf('%sdeployments/%s_SMS_mapSH.png',conf.mapsdir,EXP),...
                [folder_output_SMS NATION '/MAPS/'],'f');
            
            data1 = readtable('info_tags_sms.csv');
            data2 = data1(strcmp(data1.group,NATION),:);
            data3 = data2(data2.ispublic==1,:);
            writetable(data3,['info_tags_sms_' NATION '.csv']);
            [status,message] = movefile(['info_tags_sms_' NATION '.csv'],[folder_output_SMS NATION],'f');
            
        end
    end
end

% zip files
disp('Zip folders');    
folder_output = sprintf('%s%s/',conf.public,conf.version_SMS);
zip(sprintf('%s%s_ALL.zip',conf.public,conf.version_SMS),sprintf('%s%s',conf.public,conf.version_SMS));

