%% copy final data files in public folder

folder_output = sprintf('%s%s/',conf.public,conf.version);
[status,msg,msgID] = rmdir(folder_output, 's');
[s,mess,messid] = mkdir(folder_output);


[status,message] = copyfile('README_licenseODbl.txt',folder_output,'f');
[status,message] = copyfile([conf.processdir 'info_total.csv'],folder_output,'f');
[status,message] = copyfile([conf.processdir 'info_tags.csv'],folder_output,'f');
[status,message] = copyfile([conf.processdir 'info_deployments.csv'],folder_output,'f');
[status,message] = copyfile([conf.processdir 'info_groups.csv'],folder_output,'f');
[status,message] = copyfile(sprintf('%sglobal/map_global_public.png',conf.mapsdir),folder_output,'f');
[status,message] = copyfile(sprintf('%sglobal/map_SH.png',conf.mapsdir),folder_output,'f');
[status,message] = copyfile(sprintf('%sglobal/map_NH.png',conf.mapsdir),folder_output,'f');

EXP_all=tags_processed(conf);
list_NATION=unique(EXP_all.country);

for  kNATION=1:length(list_NATION),
    
    NATION = list_NATION{kNATION};
    EXPs=tags_processed(conf,NATION);
    list_EXP = EXPs.deployment_code;
    
    if ~any(EXPs.process==1), continue; end
    if ~any(EXPs.public ==1), continue; end
    
    [s,mess,messid] = mkdir([folder_output NATION '/PDF/']);
    [s,mess,messid] = mkdir([folder_output NATION '/MAPS/']);
    [s,mess,messid] = mkdir([folder_output NATION '/DATA_ncARGO/']);
    [s,mess,messid] = mkdir([folder_output NATION '/DATA_csv_interp/']);
    [s,mess,messid] = mkdir([folder_output NATION '/METADATA/']);
    
    data1 = readtable([conf.processdir 'info_tags.csv']);
    data2 = data1(strcmp(data1.group,NATION),:);
    data3 = data2(data2.ispublic==1,:);
    writetable(data3,[folder_output NATION '/info_tags_' NATION '.csv']);
    
    data1 = readtable([conf.processdir 'info_deployments.csv']);
    data2 = data1(strcmp(data1.group,NATION),:);
    data3 = data2(data2.ispublic==1,:);
    writetable(data3,[folder_output NATION '/info_deployments_' NATION '.csv']);
    
    [status,message] = copyfile('README_licenseODbl.txt',[folder_output NATION],'f');
    [status,message] = copyfile(sprintf('%sgroups/%s_*.png',conf.mapsdir,NATION),[folder_output NATION],'f');

    for kEXP = 1:length(list_EXP),
        
        EXP = list_EXP{kEXP};
        info_deployment=load_info_deployment(conf,EXP);
        ispublic = logical(info_deployment.public);
        if ispublic
            [status,message] = copyfile( ...
                sprintf('%s%s_lr1_doc_adj.pdf',conf.texdir,EXP),...
                [folder_output NATION '/PDF/'],'f');
            try
                [status,message] = copyfile( ...
                    sprintf('%s%s_hr2_doc_adj.pdf',conf.texdir,EXP),...
                    [folder_output NATION '/PDF/'],'f');
            end
            [status,message] = copyfile( ...
                sprintf('%sdeployments/%s_map*.png',conf.mapsdir,EXP),...
                [folder_output NATION '/MAPS/'],'f');
            
            list_tag = info_deployment.list_tag_lr1;
            Ntag = length(list_tag);
            for jj=1:Ntag,
                
                ncfile = sprintf('%s%s',info_deployment.dir,list_tag(jj).name);
                smru_name = list_tag(jj).name(1:end-12);
                ncfile_final = strrep(list_tag(jj).name, '_lr1', '');
                ncfile_hr1 = strrep(ncfile, '_lr1', '_hr1');
                ncfile_hr2 = strrep(ncfile, '_lr1', '_hr2');
                odvfile  = sprintf('%s%s_ODV.txt',info_deployment.dir,smru_name);
                odvfile2 = sprintf('%s%s_ODV.txt.zip',info_deployment.dir,smru_name);
                metafile = sprintf('%s%s_METADATA.txt',info_deployment.dir,smru_name);
                metafile2 = sprintf('%s%s_METADATA.json',info_deployment.dir,smru_name);
                
                if exist(ncfile_hr2,'file')
                    [status,message] = copyfile(ncfile_hr2,[folder_output NATION '/DATA_ncARGO/' ncfile_final],'f');
                else
                    [status,message] = copyfile(ncfile,[folder_output NATION '/DATA_ncARGO/' ncfile_final],'f');
                end
                
                if exist(odvfile,'file'),
                    zip(odvfile2,odvfile);
                    [status,message] = movefile(odvfile2,[folder_output NATION '/DATA_csv_interp/' ],'f');
                    [status,message] = copyfile(metafile,[folder_output NATION '/METADATA/' ],'f');
                    [status,message] = copyfile(metafile2,[folder_output NATION '/METADATA/' ],'f');
                end
                
            end
            
            %% update list_deployment
            if isempty(conf.list_deployment{EXP,'first_version'}{1})
                conf.list_deployment{EXP,'first_version'} = {conf.version};
            end
            conf.list_deployment{EXP,'last_version'} = {conf.version};
            writetable(conf.list_deployment,[conf.processdir 'list_deployment.csv'],'WriteRowNames',true);

        end
    end

end


% zip files
disp('Zip folders');    
folder_output = sprintf('%s%s/',conf.public,conf.version);
EXP_all=tags_processed(conf);
list_NATION=unique(EXP_all.country);
for kk = 1:length(list_NATION)
    NATION = list_NATION{kk};
    if isfolder(sprintf('%s%s/%s',conf.public,conf.version,NATION)),
        zip(sprintf('%s%s_%s.zip',conf.public,conf.version,NATION),...
            sprintf('%s%s/%s',conf.public,conf.version,NATION));
    end
end
zip(sprintf('%s%s_ALL.zip',conf.public,conf.version),sprintf('%s%s',conf.public,conf.version));

