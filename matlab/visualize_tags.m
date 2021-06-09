function visualize_tags(conf,EXP)

if isempty(conf),
    conf = init_mirounga;
end
info_deployment=load_info_deployment(conf,EXP);
if isempty(info_deployment.list_tag),
    return
end

[SUCCESS,MESSAGE,MESSAGEID] = copyfile([conf.texdir EXP '*.pdf'],[conf.plotdir EXP]);

if exist([conf.calibplotdir EXP],'dir')
    d = dir([conf.calibplotdir EXP '/calibration_*']);
    if ~isempty(d),
        for kk = 1:length(d),
            nfile = d(kk).name;
            nfile2 = strrep(nfile,'calibration_','');
            [SUCCESS,MESSAGE,MESSAGEID] = copyfile([conf.calibplotdir EXP '/' nfile],[conf.plotdir EXP '/' nfile2]);
        end
    end
end

[SUCCESS,MESSAGE,MESSAGEID] = copyfile([conf.mapsdir 'deployments/' EXP '*.png'],[conf.plotdir EXP]);
[SUCCESS,MESSAGE,MESSAGEID] = copyfile([conf.datadir EXP '/*_METADATA.txt'],[conf.plotdir EXP]);
[SUCCESS,MESSAGE,MESSAGEID] = copyfile([conf.datadir EXP '/*_diary.txt'],[conf.plotdir EXP]);

for kk=1:length(info_deployment.list_smru_name)
    smru_name = info_deployment.list_smru_name{kk};
    if ~ismember(smru_name,conf.table_filter.smru_platform_name),
        conf.table_filter(end+1,:)={smru_name,0,'date_min',0,NaN};
    end
end

name_file=[conf.processdir 'table_filter.csv'];
writetable(conf.table_filter,name_file,'WriteRowNames',1,'Delimiter',',');
