function [table_depl,table_tag] = share_meop_data(list_EXP,sharedir)
%% copy final data files in public folder

    init_config;
    conf = init_mirounga;
    if isstr(list_EXP) | iscell(list_EXP),
        list_EXP = tags_processed(list_EXP);
    end
    if isempty(list_EXP),
        return
    end
    
    [status,msg,msgID] = rmdir(sharedir, 's');
    [s,mess,messid] = mkdir(sharedir);
    if ~strcmp(sharedir(end),'/'),
        sharedir(end+1)='/';
    end

    table_depl_empty  = table({},{},[],[],[],[],...
        'VariableNames',{'depl','group','ispublic','Ntag','NprofTS','NprofT'});
    table_tag_empty   = table({},{},{},[],[],[],{},{},{},{},{},{},{},{},{},...
        'VariableNames',{'tag','depl','group','ispublic','NprofTS','NprofT',...
        'start_date','PTT','WMO','doi','body','platform_code',...
        'location','species','location_class'});

    table_depl = table_depl_empty;
    table_tag = table_tag_empty;
    for kEXP = 1:size(list_EXP,1),
    
        EXP = list_EXP.Properties.RowNames{kEXP};
        info_deployment=load_info_deployment(conf,EXP);
        ispublic = logical(info_deployment.public);        
        NATION = info_deployment.NATION;
        table_aux = cell2table({EXP,NATION,ispublic,0,0,0});
        table_aux.Properties.VariableNames = table_depl.Properties.VariableNames;
        table_depl = [ table_depl ; table_aux ];

        [s,mess,messid] = mkdir([sharedir EXP]);
        try
            [status,message] = copyfile( ...
                sprintf('%s%s_hr2_doc_adj.pdf',conf.texdir,EXP),...
                [sharedir EXP],'f');
        catch
            [status,message] = copyfile( ...
                sprintf('%s%s_lr1_doc_adj.pdf',conf.texdir,EXP),...
                [sharedir EXP],'f');

        end
        [status,message] = copyfile( ...
            sprintf('%sdeployments/%s_map*.png',conf.mapsdir,EXP),...
            [sharedir EXP],'f');
        
        list_tag = info_deployment.list_tag_hr1;
        table_tag_aux = table_tag_empty;
        for jj=1:length(list_tag),
            
            name_prof = [info_deployment.dir list_tag(jj).name];
            [smru_prefix,Nsplit,suffix,namedir] = smru_name_from_name_prof(name_prof);
            smru_name = gen_smru_name(smru_prefix,Nsplit);
            ncfile_hr1 = gen_name_prof(smru_prefix,Nsplit,'hr1',namedir);
            ncfile_hr2 = gen_name_prof(smru_prefix,Nsplit,'hr2',namedir);
            ncfile_lr1 = gen_name_prof(smru_prefix,Nsplit,'lr1',namedir);

            ncfile_final = gen_name_prof(smru_prefix,Nsplit,'',[sharedir EXP]);
            ncfile_raw = gen_name_prof(smru_prefix,Nsplit,'raw',[sharedir EXP]);
            metafile = sprintf('%s%s_METADATA.txt',info_deployment.dir,smru_name);
            metafile2 = sprintf('%s%s_METADATA.json',info_deployment.dir,smru_name);
            
            M=ARGO_load_qc(name_prof,1);
            Mattr=ncloadatt_struct(name_prof);
            M.Tmask=double(M.TEMP_QC<2);
            M.Smask=double(M.PSAL_QC<2);
            
            % update tables
            NprofTS=length(find(sum(M.Tmask.*M.Smask)~=0));
            NprofT =length(find(sum(M.Tmask)~=0));
            
            table_depl{end,4:6}   = table_depl{end,4:6} + [1 NprofTS NprofT];
            year = datestr(min(M.JULD(find(sum(M.Tmask)'))),29);
            if ~isfield(Mattr,'loc_algorithm'), Mattr.loc_algorithm = Mattr.location_class; end
            table_aux = cell2table({M.smru_platform_code,EXP,NATION,ispublic,NprofTS,NprofT...
                ,year,Mattr.ptt,Mattr.wmo_platform_code,Mattr.reference_doi,...
                Mattr.instr_id,Mattr.platform_code,Mattr.location,...
                Mattr.species,Mattr.loc_algorithm});
            table_aux.Properties.VariableNames = table_tag.Properties.VariableNames;
            table_tag_aux = [ table_tag_aux ; table_aux ];

        
            if exist(ncfile_hr2,'file')
                [status,message] = copyfile(ncfile_hr2,ncfile_final,'f');
            else
                [status,message] = copyfile(ncfile_hr1,ncfile_final,'f');
            end
            [status,message] = copyfile(ncfile_lr1,ncfile_raw,'f');
            [status,message] = copyfile(metafile2,[sharedir EXP],'f');    

        end
        writetable(table_tag_aux   ,[sharedir EXP '/info_tags.csv']);
        table_tag = [ table_tag ; table_tag_aux ];
    end

    writetable(table_depl  ,[sharedir 'info_deployments.csv']);
    writetable(table_tag   ,[sharedir 'info_tags.csv']);
    [status,message] = copyfile('README_licenseODbl.txt',sharedir,'f');






