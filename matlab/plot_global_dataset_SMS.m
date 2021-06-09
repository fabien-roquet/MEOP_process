function plot_global_dataset_SMS(conf)
%% perform plots and pdf document on the general MEOP-CTD dataset

if isempty(conf),
    conf = init_mirounga;
end

do_plot_figure = 1;

close all

% sort deployments by NATION
EXP_all=tags_processed(conf);
list_NATION=unique(EXP_all.country);

%% plot map all

LAT=[]; LON=[]; DATE=[];
Nprof1=0;Ntag1=0;H1=[];nation1=[];Ndep1=[];
Nprof2=0;Ntag2=0;H2=[];nation2=[];Ndep2=[];
Nprof3=0;Ntag3=0;H3=[];nation3=[];Ndep3=[];
Nprof4=0;Ntag4=0;H4=[];nation4=[];Ndep4=[];name_tag4 = {};
Nprof5=0;Ntag5=0;H5=[];nation5=[];Ndep5=[];name_tag5 = {};
Nprof6=0;Ntag6=0;H6=[];nation6=[];Ndep6=[];
Nprof7=0;Ntag7=0;H7=[];nation7=[];Ndep7=[];
Nprof8=0;Ntag8=0;H8=[];nation8=[];Ndep8=[];
Nprof9=0;Ntag9=0;H9=[];nation9=[];Ndep9=[];
Nprof10=0;Ntag10=0;H10=[];nation10=[];Ndep10=[];
Nprof11=0;Ntag11=0;H11=[];nation11=[];Ndep11=[];
Nprof12=0;Ntag12=0;H12=[];nation12=[];Ndep12=[];name_tag12 = {};

table_tag_hr= table({},{},{},[],[],[],{},{},{},{},{},{},{},{},{},'VariableNames',{'tag','depl','group','ispublic','NprofTS','NprofT','start_date','PTT','WMO','doi','body','platform_code','location','species','location_class'});

if do_plot_figure
    
    set(0, 'DefaultFigureVisible', 'off');

    % group SH, continuous data
    hf(11) = figure(11);clf,
    set(hf(11), 'Color', 'w');
    set(hf(11), 'units','centimeters','Position', [25 5 20 20]);
    m_proj('stereographic','lat',-90,'long',90,'radius',60);
    m_grid('box','on','xtick',12,...
        'tickdir','out','ytick',[-90:10:-40],...
        'XaxisLocation','top','linest','-');
    topoplot([-180 180 -90 -30],'m_map')
    hold on
    % deployment SH, continuous data
    hf(12) = figure(12);clf,
    set(hf(12), 'Color', 'w');
    set(hf(12), 'units','centimeters','Position', [25 5 20 20]);
    m_proj('stereographic','lat',-90,'long',90,'radius',60);
    m_grid('box','on','xtick',12,...
        'tickdir','out','ytick',[-90:10:-40],...
        'XaxisLocation','top','linest','-');
    topoplot([-180 180 -90 -30],'m_map')
    hold on
    %
end

for  kNATION=1:length(list_NATION),
    
    NATION = list_NATION{kNATION};
    EXPs=tags_processed(conf,NATION);
    list_EXP = EXPs.deployment_code;
    delete(H11); Nprof11=0;Ntag11=0;H11=[];nation11=[];Ndep11=[];
    col_group=linspecer(length(list_EXP)); numcolgroup=0;
    if ~any(EXPs.process==1), continue; end
    
    for kEXP = 1:length(list_EXP),
        
        EXP = list_EXP{kEXP};
        numcolgroup=numcolgroup+1;
        info_deployment=load_info_deployment(conf,EXP);
        list_tag = info_deployment.list_tag_fr1;
        Ntag = length(list_tag);
        
        name_tag12 = {}; delete(H12); Nprof12=0;Ntag12=0;H12=[];nation12=[];Ndep12=[];
        public=logical(info_deployment.public);
        
        for jj=1:Ntag,
            
            name_prof = sprintf('%s%s',info_deployment.dir,list_tag(jj).name);
            smru_name = list_tag(jj).name(1:end-12);
            exist_continuous = conf.hr_continuous(find(strcmp(smru_name,conf.hr_smru_name)));
            
            if ~exist_continuous, continue; end
            
            ncfile_fr1 = strrep(name_prof, '_lr0', '_fr1');
            Mhr=ARGO_load_qc(ncfile_fr1,1);
            Mattr=ncloadatt_struct(name_prof);
            Mhr.Tmask=double(Mhr.TEMP_QC<2);
            Mhr.Smask=double(Mhr.PSAL_QC<2);
            NprofTS_hr=length(find(sum(Mhr.Tmask.*Mhr.Smask)~=0));
            NprofT_hr =length(find(sum(Mhr.Tmask)~=0));
            
            if NprofT_hr==0, continue, end
            
            year = datestr(min(Mhr.JULD(find(sum(Mhr.Tmask)'))),29);
            table_aux = cell2table({Mhr.smru_platform_code,EXP,NATION,public,NprofTS_hr,NprofT_hr...
                ,year,Mattr.ptt,Mattr.wmo_platform_code,Mattr.reference_doi,Mattr.instr_id,Mattr.platform_code,Mattr.location,Mattr.species,Mattr.loc_algorithm});
            table_aux.Properties.VariableNames = table_tag_hr.Properties.VariableNames;
            table_tag_hr = [ table_tag_hr ; table_aux ];
            
            if Mhr.ntag~=1, error('Should be only one tag per file'); end
            descr=Mhr.list_descr{1};
            Khr=find(strcmp(Mhr.platform_number,descr) & sum(Mhr.Tmask)'>0);
            if length(Khr)==0, continue, end

            disp(['process tag: ' smru_name]);
            
            if do_plot_figure
                
                % 11 & 12
                if nanmean(Mhr.LATITUDE)<-20
                    
                    % fig 11
                    Ntag11=Ntag11+1;
                    Nprof11=Nprof11+length(Khr);
                    nation11{Ntag11,1}=info_deployment.NATION;
                    Ndep11{Ntag11,1}=info_deployment.EXP;
                    figure(11); m_proj('stereographic','lat',-90,'long',90,'radius',60); hold on
                    h=m_plot(Mhr.LONGITUDE(Khr),Mhr.LATITUDE(Khr));
                    col= col_group(numcolgroup,:);
                    set(h,'linewidth',.5,'color',col);
                    H11=[H11 h(1)];
                    % fig 12
                    Ntag12=Ntag12+1;
                    Nprof12=Nprof12+length(Khr);
                    nation12{Ntag12,1}=info_deployment.NATION;
                    Ndep12{Ntag12,1}=info_deployment.EXP;
                    name_tag12{Ntag12} = Mhr.smru_platform_code;
                    figure(12); m_proj('stereographic','lat',-90,'long',90,'radius',60); hold on
                    h=m_plot(Mhr.LONGITUDE(Khr),Mhr.LATITUDE(Khr));
                    set(h,'linewidth',.5);
                    H12=[H12 h(1)];
                    
                end
                
            end
            
        end
        
        if do_plot_figure
            if Ntag12,
                [NATION12,Icol12]=unique(nation12,'stable');
                for kk=1:length(name_tag12),
                    name_tag12{kk}=strrep(name_tag12{kk},'_','\_');
                end
                set(0, 'currentfigure', hf(12));
                suptitle(sprintf('Deployment %s (%s), SMS data : %d profiles, %d tags',info_deployment.EXP,NATION12{1},Nprof12,Ntag12));
                h=legend(H12,name_tag12,'fontsize',8,'Location','bestoutside'); H12=[H12 h];
                nfile=sprintf('%sdeployments/%s_SMS_mapSH',conf.mapsdir,strtrim(EXP));
                print(nfile,'-dpng','-r300');
                delete(H12); H12=[];
            end
        end
        
    end
    
    if do_plot_figure
        if Ntag11,
            [name_dep11,Icol11]=unique(Ndep11,'stable');
            set(0, 'currentfigure', hf(11));
            suptitle(sprintf('Group %s : %d profiles, %d deployments, %d tags',NATION,Nprof11,length(name_dep11),Ntag11));
            h=legend(H11(Icol11),name_dep11,'fontsize',8,'Location','bestoutside'); H11=[H11 h];
            nfile=sprintf('%sgroups/%s_SMS_mapSH',conf.mapsdir,strtrim(NATION));
            print(nfile,'-dpng','-r300');
            delete(H11); H11=[];
        end
    end
    
end

% save tables
writetable(table_tag_hr,[conf.processdir 'info_tags_sms.csv']);


