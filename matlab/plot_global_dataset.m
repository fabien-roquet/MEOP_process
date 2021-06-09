function plot_global_dataset(conf)

if isempty(conf),
    conf = init_mirounga;
end


%% perform plots and pdf document on the general MEOP-CTD dataset

do_plot_figure = 1;
do_density_plot = 0;

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

table_total = table({'total';'public';'private'},[0;0;0],[0;0;0],[0;0;0],[0;0;0],[0;0;0],'VariableNames',{'type','Ngroup','Ndepl','Ntag','NprofTS','NprofT'});
table_group = table({},[],[],[],[],[],[],[],[],'VariableNames',{'group','Ndepl_public','Ntag_public','NprofTS_public','NprofT_public','Ndepl_private','Ntag_private','NprofTS_private','NprofT_private'});
table_depl  = table({},{},[],[],[],[],'VariableNames',{'depl','group','ispublic','Ntag','NprofTS','NprofT'});
table_tag   = table({},{},{},[],[],[],{},{},{},{},{},{},{},{},{},'VariableNames',{'tag','depl','group','ispublic','NprofTS','NprofT','start_date','PTT','WMO','doi','body','platform_code','location','species','location_class'});

if do_plot_figure

    set(0, 'DefaultFigureVisible', 'off');
    
    % total SH, ispublic
    fh(1) = figure(1);clf,
    set(fh(1), 'Color', 'w');
    set(fh(1), 'units','centimeters','Position', [2 5 20 20]);
    m_proj('stereographic','lat',-90,'long',90,'radius',60);
    m_grid('box','on','xtick',12,...
        'tickdir','out','ytick',[-90:10:-40],...
        'XaxisLocation','top','linest','-');
    topoplot([-180 180 -90 -30],'m_map')
    hold on
    % total NH, ispublic
    fh(2) = figure(2);clf,
    set(fh(2), 'Color', 'w');
    set(fh(2), 'units','centimeters','Position', [15 15 20 20]);
    m_proj('stereographic','lat',90,'long',90,'radius',60);
    m_grid('box','on','xtick',12,...
        'tickdir','out','ytick',[40:10:90],...
        'XaxisLocation','bottom','linest','-');
    topoplot([-180 180 30 90],'m_map')
    hold on
    % total by nation, ispublic
    fh(3) = figure(3);clf,
    set(fh(3), 'Color', 'w');
    set(fh(3), 'units','centimeters','Position', [5 30 20 14]);
    topoplot([-250 110 -80 80]); latratio(45)
    set(gca,'xtick',[-240:60:60],'xticklabel',{'120E','180E','120W','60W','0E','60E'});
    set(gca,'ytick',[-80:20:80],'yticklabel',{'80S','60S','40S','20S','0N','20N','40N','60N','80N'});
    grid on,    hold on
    % deployment SH
    fh(4) = figure(4);clf,
    set(fh(4), 'Color', 'w');
    set(fh(4), 'units','centimeters','Position', [25 5 20 20]);
    m_proj('stereographic','lat',-90,'long',90,'radius',60);
    m_grid('box','on','xtick',12,...
        'tickdir','out','ytick',[-90:10:-40],...
        'XaxisLocation','top','linest','-');
    topoplot([-180 180 -90 -30],'m_map')
    hold on
    % deployment NH
    fh(5) = figure(5);clf,
    set(fh(5), 'Color', 'w');
    set(fh(5), 'units','centimeters','Position', [35 25 20 20]);
    m_proj('stereographic','lat',90,'long',90,'radius',60);
    m_grid('box','on','xtick',12,...
        'tickdir','out','ytick',[40:10:90],...
        'XaxisLocation','bottom','linest','-');
    topoplot([-180 180 30 90],'m_map')
    hold on
    % group SH, ispublic
    fh(6) = figure(6);clf,
    set(fh(6), 'Color', 'w');
    set(fh(6), 'units','centimeters','Position', [50 5 20 20]);
    m_proj('stereographic','lat',-90,'long',90,'radius',60);
    m_grid('box','on','xtick',12,...
        'tickdir','out','ytick',[-90:10:-40],...
        'XaxisLocation','top','linest','-');
    topoplot([-180 180 -90 -30],'m_map')
    hold on
    % group NH, ispublic
    fh(7) = figure(7);clf,
    set(fh(7), 'Color', 'w');
    set(fh(7), 'units','centimeters','Position', [60 25 20 20]);
    m_proj('stereographic','lat',90,'long',90,'radius',60);
    m_grid('box','on','xtick',12,...
        'tickdir','out','ytick',[40:10:90],...
        'XaxisLocation','bottom','linest','-');
    topoplot([-180 180 30 90],'m_map')
    hold on
    % total public
    fh(8) = figure(8);clf,
    set(fh(8), 'Color', 'w');
    set(fh(8), 'units','centimeters','Position', [5 30 20 14]);
    topoplot([-250 110 -80 80]); latratio(45)
    set(gca,'xtick',[-240:60:60],'xticklabel',{'120E','180E','120W','60W','0E','60E'});
    set(gca,'ytick',[-80:20:80],'yticklabel',{'80S','60S','40S','20S','0N','20N','40N','60N','80N'});
    grid on,    hold on
    % total private
    fh(9) = figure(9);clf,
    set(fh(9), 'Color', 'w');
    set(fh(9), 'units','centimeters','Position', [5 30 20 14]);
    topoplot([-250 110 -80 80]); latratio(45)
    set(gca,'xtick',[-240:60:60],'xticklabel',{'120E','180E','120W','60W','0E','60E'});
    set(gca,'ytick',[-80:20:80],'yticklabel',{'80S','60S','40S','20S','0N','20N','40N','60N','80N'});
    grid on,    hold on
    % total SH ellies-only
    fh(10) = figure(10);clf,
    set(fh(10), 'Color', 'w');
    set(fh(10), 'units','centimeters','Position', [2 5 20 20]);
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
    if ~any(EXPs.process==1), continue; end

    delete(H6); Nprof6=0;Ntag6=0;H6=[];nation6=[];Ndep6=[];
    delete(H7); Nprof7=0;Ntag7=0;H7=[];nation7=[];Ndep7=[];
    col_group=linspecer(length(list_EXP)); numcolgroup=0;
    table_total{1:2,2} = table_total{1:2,2} + 1;
    if any(EXPs.public==1), table_total{3,2} = table_total{3,2} + 1; end
    table_aux = cell2table({NATION,0,0,0,0,0,0,0,0});
    table_aux.Properties.VariableNames = table_group.Properties.VariableNames;
    table_group = [ table_group ; table_aux ];
    
    for kEXP = 1:length(list_EXP),
        
        EXP = list_EXP{kEXP};
        numcolgroup=numcolgroup+1;
        disp(['process deployment: ' EXP]);
        info_deployment=load_info_deployment(conf,EXP);
        list_tag = dir([info_deployment.dir '*_hr1_prof.nc']);
        Ntag = length(list_tag);
        name_tag4 = {}; delete(H4); Nprof4=0;Ntag4=0;H4=[];nation4=[];Ndep4=[];
        name_tag5 = {}; delete(H5); Nprof5=0;Ntag5=0;H5=[];nation5=[];Ndep5=[];
        ispublic=logical(info_deployment.public);
        table_total{1,3} = table_total{1,3} + 1;
        if ispublic,
            table_total{2,3} = table_total{2,3} + 1;
            table_group{end,2} = table_group{end,2} + 1;
        else,
            table_total{3,3} = table_total{3,3} + 1;
            table_group{end,6} = table_group{end,6} + 1;
        end
        table_aux = cell2table({EXP,NATION,ispublic,0,0,0});
        table_aux.Properties.VariableNames = table_depl.Properties.VariableNames;
        table_depl = [ table_depl ; table_aux ];
        
        for jj=1:Ntag,
            
            name_prof = sprintf('%s%s',info_deployment.dir,list_tag(jj).name);
            smru_name = list_tag(jj).name(1:end-12);
            
            M=ARGO_load_qc(name_prof,1);
            Mattr=ncloadatt_struct(name_prof);
            M.Tmask=double(M.TEMP_QC<2);
            M.Smask=double(M.PSAL_QC<2);
            
            % update tables
            NprofTS=length(find(sum(M.Tmask.*M.Smask)~=0));
            NprofT =length(find(sum(M.Tmask)~=0));
            
            if NprofT==0, continue, end
            
            table_total{1,4:6} = table_total{1,4:6} + [1 NprofTS NprofT];
            if ispublic,
                table_total{2,4:6}   = table_total{2,4:6} + [1 NprofTS NprofT];
                table_group{end,3:5} = table_group{end,3:5} + [1 NprofTS NprofT];
            else
                table_total{3,4:6} = table_total{3,4:6} + [1 NprofTS NprofT];
                table_group{end,7:9} = table_group{end,7:9} + [1 NprofTS NprofT];
            end
            table_depl{end,4:6}   = table_depl{end,4:6} + [1 NprofTS NprofT];
            year = datestr(min(M.JULD(find(sum(M.Tmask)'))),29);
            if ~isfield(Mattr,'loc_algorithm'), Mattr.loc_algorithm = Mattr.location_class; end
            table_aux = cell2table({M.smru_platform_code,EXP,NATION,ispublic,NprofTS,NprofT...
                ,year,Mattr.ptt,Mattr.wmo_platform_code,Mattr.reference_doi,...
                Mattr.instr_id,Mattr.platform_code,Mattr.location,...
                Mattr.species,Mattr.loc_algorithm});
            table_aux.Properties.VariableNames = table_tag.Properties.VariableNames;
            table_tag = [ table_tag ; table_aux ];
            
            if M.ntag~=1, error('Should be only one tag per file'); end
            descr=M.list_descr{1};
            K=find(strcmp(M.platform_number,descr) & sum(M.Tmask)'>0);
            
            lon=M.LONGITUDE(K); lat=M.LATITUDE(K);
            lon(lon>180)=lon(lon>180)-360;
            LAT=[LAT;lat]; LON=[LON;lon]; DATE=[DATE;M.JULD_LOCATION(K)];
            
            if do_plot_figure
                
                if length(K)

                    % fig 1 & 4 & 6 & 10 & 11 & 12
                    if nanmean(M.LATITUDE)<-20
                    
                        % fig 1
                        if ispublic
                            Ntag1=Ntag1+1;
                            Nprof1=Nprof1+length(K);
                            nation1{Ntag1,1}=info_deployment.NATION;
                            Ndep1{Ntag1,1}=info_deployment.EXP;
                            set(0, 'currentfigure', fh(1)); m_proj('stereographic','lat',-90,'long',90,'radius',60); hold on
                            h=m_plot(M.LONGITUDE(K),M.LATITUDE(K));
                            col= conf.list_color(find(strcmp(NATION,conf.list_group)),:);
                            set(h,'linewidth',.3,'color',col);
                            H1=[H1 h(1)];
                        end
                        % fig 4
                        Ntag4=Ntag4+1;
                        Nprof4=Nprof4+length(K);
                        nation4{Ntag4,1}=info_deployment.NATION;
                        Ndep4{Ntag4,1}=info_deployment.EXP;
                        name_tag4{Ntag4} = M.smru_platform_code;
                        set(0, 'currentfigure', fh(4)); m_proj('stereographic','lat',-90,'long',90,'radius',60); hold on
                        h=m_plot(M.LONGITUDE(K),M.LATITUDE(K));
                        set(h,'linewidth',.3);
                        H4=[H4 h];
                        % fig 6
                        if ispublic
                            Ntag6=Ntag6+1;
                            Nprof6=Nprof6+length(K);
                            nation6{Ntag6,1}=info_deployment.NATION;
                            Ndep6{Ntag6,1}=info_deployment.EXP;
                            set(0, 'currentfigure', fh(6)); m_proj('stereographic','lat',-90,'long',90,'radius',60); hold on
                            h=m_plot(M.LONGITUDE(K),M.LATITUDE(K));
                            col= col_group(numcolgroup,:);
                            set(h,'linewidth',.3,'color',col);
                            H6=[H6 h];
                        end
                        % fig 10
                        if ismember(Mattr.species,{'Southern ellie','Ellies'})
                            Ntag10=Ntag10+1;
                            Nprof10=Nprof10+length(K);
                            nation10{Ntag10,1}=info_deployment.NATION;
                            Ndep10{Ntag10,1}=info_deployment.EXP;
                            set(0, 'currentfigure', fh(10)); m_proj('stereographic','lat',-90,'long',90,'radius',60); hold on
                            h=m_plot(M.LONGITUDE(K),M.LATITUDE(K));
                            col= conf.list_color(find(strcmp(NATION,conf.list_group)),:);
                            set(h,'linewidth',.3,'color',col);
                            H10=[H10 h(1)];
                        end
                    end
                    
                    % fig 2 & 5 & 7
                    if nanmean(M.LATITUDE)>20
                        % fig 2
                        if ispublic
                            Ntag2=Ntag2+1;
                            Nprof2=Nprof2+length(K);
                            nation2{Ntag2,1}=info_deployment.NATION;
                            Ndep2{Ntag2,1}=info_deployment.EXP;
                            set(0, 'currentfigure', fh(2)); m_proj('stereographic','lat',90,'long',90,'radius',60); hold on
                            h=m_plot(M.LONGITUDE(K),M.LATITUDE(K));
                            col= conf.list_color(find(strcmp(NATION,conf.list_group)),:);
                            set(h,'linewidth',.3,'color',col);
                                H2=[H2 h(1)];
                            end
                            % fig 5
                            Ntag5=Ntag5+1;
                            Nprof5=Nprof5+length(K);
                            nation5{Ntag5,1}=info_deployment.NATION;
                            Ndep5{Ntag5,1}=info_deployment.EXP;
                            name_tag5{Ntag5} = M.smru_platform_code;
                            set(0, 'currentfigure', fh(5)); m_proj('stereographic','lat',90,'long',90,'radius',60); hold on
                            h=m_plot(M.LONGITUDE(K),M.LATITUDE(K));
                            set(h,'linewidth',.3);
                            H5=[H5 h];
                            % fig 7
                            if ispublic
                                Ntag7=Ntag7+1;
                                Nprof7=Nprof7+length(K);
                                nation7{Ntag7,1}=info_deployment.NATION;
                                Ndep7{Ntag7,1}=info_deployment.EXP;
                                set(0, 'currentfigure', fh(7)); m_proj('stereographic','lat',90,'long',90,'radius',60); hold on
                                h=m_plot(M.LONGITUDE(K),M.LATITUDE(K));
                                col= col_group(numcolgroup,:);
                                set(h,'linewidth',.3,'color',col);
                                H7=[H7 h];
                            end
                    end
                    
                    % fig 3 & 8 & 9 & 10
                    if ispublic
                        
                        % fig 3
                        Ntag3=Ntag3+1;
                        Nprof3=Nprof3+length(K);
                        nation3{Ntag3,1}=info_deployment.NATION;
                        Ndep3{Ntag3,1}=info_deployment.EXP;
                        set(0, 'currentfigure', fh(3));
                        col= conf.list_color(find(strcmp(NATION,conf.list_group)),:);
                        h1=[]; h2=[];
                        lat2=lat; lon2=lon; lon2(lon2>110)=lon2(lon2>110)-360;
                        if any(lon2<-160) & any(lon2>20),
                            h1=plot(lon2(lon2<-70),lat2(lon2<-70)); set(h1,'linewidth',.3,'color',col);
                            h2=plot(lon2(lon2>-70),lat2(lon2>-70)); set(h2,'linewidth',.3,'color',col);
                        else
                            h1=plot(lon2,lat2); set(h1,'linewidth',.3,'color',col);
                        end
                        if ~isempty(h1)
                            H3=[H3 h1(1)];
                        elseif ~isempty(h2)
                            H3=[H3 h2(1)];
                        end
                        
                        % fig 8
                        Ntag8=Ntag8+1;
                        Nprof8=Nprof8+length(K);
                        nation8{Ntag8,1}=info_deployment.NATION;
                        Ndep8{Ntag8,1}=info_deployment.EXP;
                        set(0, 'currentfigure', fh(8));
                        col= conf.list_color(find(strcmp(NATION,conf.list_group)),:);
                        h1=[]; h2=[];
                        lat2=lat; lon2=lon; lon2(lon2>110)=lon2(lon2>110)-360;
                        if any(lon2<-160) & any(lon2>20),
                            h1=plot(lon2(lon2<-70),lat2(lon2<-70)); set(h1,'linewidth',.3,'color',col);
                            h2=plot(lon2(lon2>-70),lat2(lon2>-70)); set(h2,'linewidth',.3,'color',col);
                        else
                            h1=plot(lon2,lat2); set(h1,'linewidth',.3,'color',col);
                        end
                        if ~isempty(h1)
                            H8=[H8 h1(1)];
                        elseif ~isempty(h2)
                            H8=[H8 h2(1)];
                        end
                        
                    else
                        
                        % fig 9
                        Ntag9=Ntag9+1;
                        Nprof9=Nprof9+length(K);
                        nation9{Ntag9,1}=info_deployment.NATION;
                        Ndep9{Ntag9,1}=info_deployment.EXP;
                        set(0, 'currentfigure', fh(9));
                        col= conf.list_color(find(strcmp(NATION,conf.list_group)),:);
                        h1=[]; h2=[];
                        lat2=lat; lon2=lon; lon2(lon2>110)=lon2(lon2>110)-360;
                        if any(lon2<-160) & any(lon2>20),
                            h1=plot(lon2(lon2<-70),lat2(lon2<-70)); set(h1,'linewidth',.3,'color',col);
                            h2=plot(lon2(lon2>-70),lat2(lon2>-70)); set(h2,'linewidth',.3,'color',col);
                        else
                            h1=plot(lon2,lat2); set(h1,'linewidth',.3,'color',col);
                        end
                        if ~isempty(h1)
                            H9=[H9 h1(1)];
                        elseif ~isempty(h2)
                            H9=[H9 h2(1)];
                        end
                    end
                    
                end
                
            end
            
        end
        
        if do_plot_figure
            
            if Ntag4,
                [NATION4,Icol4]=unique(nation4,'stable');
                for kk=1:length(name_tag4),
                    name_tag4{kk}=strrep(name_tag4{kk},'_','\_');
                end
                set(0, 'currentfigure', fh(4));
                suptitle(sprintf('Deployment %s (%s) : %d profiles, %d tags',info_deployment.EXP,NATION4{1},Nprof4,Ntag4));
                h=legend(H4,name_tag4,'fontsize',8,'Location','bestoutside'); H4=[H4 h];
                nfile=sprintf('%sdeployments/%s_mapSH',conf.mapsdir,strtrim(EXP));
                print(nfile,'-dpng','-r300');
                %eval(['export_fig ' nfile '.png -m5 -nocrop']);
                delete(H4); H4=[];
            end
            
            if Ntag5,
                [NATION5,Icol5]=unique(nation5,'stable');
                for kk=1:length(name_tag5),
                    name_tag5{kk}=strrep(name_tag5{kk},'_','\_');
                end
                set(0, 'currentfigure', fh(5));
                suptitle(sprintf('Deployment %s (%s) : %d profiles, %d tags',info_deployment.EXP,NATION5{1},Nprof5,Ntag5));
                h=legend(H5,name_tag5,'fontsize',8,'Location','bestoutside'); H5=[H5 h];
                nfile=sprintf('%sdeployments/%s_mapNH',conf.mapsdir,strtrim(EXP));
                print(nfile,'-dpng','-r300');
                %eval(['export_fig ' nfile '.png -m5 -nocrop']);
                delete(H5); H5=[];
            end
            
        end
        
    end
    
    if do_plot_figure
        
        if Ntag6,
            [name_dep6,Icol6]=unique(Ndep6,'stable');
            set(0, 'currentfigure', fh(6));
            suptitle(sprintf('Group %s : %d profiles, %d deployments, %d tags',NATION,Nprof6,length(name_dep6),Ntag6));
            h=legend(H6(Icol6),name_dep6,'fontsize',8,'Location','bestoutside'); H6=[H6 h];
            nfile=sprintf('%sgroups/%s_mapSH',conf.mapsdir,strtrim(NATION));
            print(nfile,'-dpng','-r300');
            %eval(['export_fig ' nfile '.png -m5 -nocrop']);
            delete(H6); H6=[];
        end
        
        if Ntag7,
            [name_dep7,Icol7]=unique(Ndep7,'stable');
            set(0, 'currentfigure', fh(7));
            suptitle(sprintf('Group %s : %d profiles, %d deployments, %d tags',NATION,Nprof7,length(name_dep7),Ntag7));
            h=legend(H7(Icol7),name_dep7,'fontsize',8,'Location','bestoutside'); H7=[H7 h];
            nfile=sprintf('%sgroups/%s_mapNH',conf.mapsdir,strtrim(NATION));
            print(nfile,'-dpng','-r300');
            %eval(['export_fig ' nfile '.png -m5 -nocrop']);
            delete(H7); H7=[];
        end
        
    end
    
end

if do_plot_figure
    
    % fig 1
    [NATION1,Icol1]=unique(nation1,'stable');
    dep1=length(unique(Ndep1));
    H1=H1(Icol1);
    set(0, 'currentfigure', fh(1)); suptitle(sprintf('MEOP-CTD SH dataset : %d profiles, %d deployments, %d tags',Nprof1,dep1,Ntag1));
    legend(H1,NATION1','fontsize',8,'Location','bestoutside');
    nfile=sprintf('%sglobal/map_SH',conf.mapsdir);
    print(nfile,'-dpng','-r300');
    %eval(['export_fig ' nfile '.png -m5 -nocrop']);
    %eval(['export_fig ' nfile '_1.pdf -m5 -nocrop -painters']);
    
    % fig 2
    [NATION2,Icol2]=unique(nation2,'stable');
    dep2=length(unique(Ndep2));
    H2=H2(Icol2);
    set(0, 'currentfigure', fh(2)); suptitle(sprintf('MEOP-CTD NH dataset : %d profiles, %d deployments, %d tags',Nprof2,dep2,Ntag2));
    legend(H2,NATION2','fontsize',8,'Location','bestoutside');
    nfile=sprintf('%sglobal/map_NH',conf.mapsdir);
    print(nfile,'-dpng','-r300');
    %eval(['export_fig ' nfile '.png -m5 -nocrop']);
    
    % fig 3
    [NATION3,Icol3]=unique(nation3,'stable');
    dep3=length(unique(Ndep3));
    H3=H3(Icol3);
    set(0, 'currentfigure', fh(3)); suptitle(sprintf('MEOP-CTD dataset : %d profiles, %d deployments, %d tags',Nprof3,dep3,Ntag3));
    legend(H3,NATION3','fontsize',8,'Location','NorthEast');
    nfile=sprintf('%sglobal/map_global',conf.mapsdir);
    print(nfile,'-dpng','-r300');
    %eval(['export_fig ' nfile '.png -m5 -nocrop']);
    
    % fig 8
    [NATION8,Icol8]=unique(nation8,'stable');
    dep8=length(unique(Ndep8));
    H8=H8(Icol8);
    set(0, 'currentfigure', fh(8)); suptitle(sprintf('MEOP-CTD public dataset : %d profiles, %d deployments, %d tags',Nprof8,dep8,Ntag8));
    legend(H8,NATION8','fontsize',8,'Location','NorthEast');
    nfile=sprintf('%sglobal/map_global_public',conf.mapsdir);
    print(nfile,'-dpng','-r300');
    %eval(['export_fig ' nfile '.png -m5 -nocrop']);
    
    % fig 9
    [NATION9,Icol9]=unique(nation9,'stable');
    dep9=length(unique(Ndep9));
    H9=H9(Icol9);
    set(0, 'currentfigure', fh(9)); suptitle(sprintf('MEOP-CTD private dataset : %d profiles, %d deployments, %d tags',Nprof9,dep9,Ntag9));
    legend(H9,NATION9','fontsize',8,'Location','NorthEast');
    nfile=sprintf('%sglobal/map_global_private',conf.mapsdir);
    print(nfile,'-dpng','-r300');
    %eval(['export_fig ' nfile '.png -m5 -nocrop']);
    
    % fig 10
    [NATION10,Icol10]=unique(nation10,'stable');
    dep10=length(unique(Ndep10));
    H10=H10(Icol10);
    set(0, 'currentfigure', fh(10)); suptitle(sprintf('MEOP-CTD SH ellies : %d profiles, %d deployments, %d tags',Nprof10,dep10,Ntag10));
    legend(H10,NATION10','fontsize',8,'Location','bestoutside');
    nfile=sprintf('%sglobal/map_SH_ellies',conf.mapsdir);
    print(nfile,'-dpng','-r300');
    %eval(['export_fig ' nfile '.png -m5 -nocrop']);
    %eval(['export_fig ' nfile '_1.pdf -m5 -nocrop -painters']);
    

end

% save tables
writetable(table_total ,[conf.processdir 'info_total.csv']);
writetable(table_group ,[conf.processdir 'info_groups.csv']);
writetable(table_depl  ,[conf.processdir 'info_deployments.csv']);
writetable(table_tag   ,[conf.processdir 'info_tags.csv']);

% plot global latex file
% sc_build_latex_global;

if do_density_plot,
    
    % %% data density distribution
    xbin=(-250:110)';
    ybin=(-90:90)';
    [XBIN,YBIN]=meshgrid(xbin(1:end-1)+.5,ybin(1:end-1)+.5);
    LAT2=LAT; LON2=LON; LON2(LON2>110)=LON2(LON2>110)-360;
    COUNT=twodhist(LON2,LAT2,xbin,ybin); COUNT(COUNT==0)=NaN;
    
    wod_ctd_path=[conf.woddir '../data_from_CTD_Collection.nc'];
    data_ctd_wod = ncload_struct(wod_ctd_path,'latitude','longitude');
    LAT_WOD=[data_ctd_wod.latitude];
    LON_WOD=[data_ctd_wod.longitude];
    LAT_WOD2=LAT_WOD; LON_WOD2=LON_WOD; LON_WOD2(LON_WOD2>110)=LON_WOD2(LON_WOD2>110)-360;
    COUNT_CTD_WOD=twodhist(LON_WOD2,LAT_WOD2,xbin,ybin); COUNT_CTD_WOD(COUNT_CTD_WOD==0)=NaN;
    
    wod_pfl_path=[conf.woddir '../data_from_PFL_Collection.nc'];
    data_pfl_wod = ncload_struct(wod_pfl_path,'latitude','longitude');
    LAT_WOD=[data_pfl_wod.latitude];
    LON_WOD=[data_pfl_wod.longitude];
    LAT_WOD2=LAT_WOD; LON_WOD2=LON_WOD; LON_WOD2(LON_WOD2>110)=LON_WOD2(LON_WOD2>110)-360;
    COUNT_PFL_WOD=twodhist(LON_WOD2,LAT_WOD2,xbin,ybin); COUNT_PFL_WOD(COUNT_PFL_WOD==0)=NaN;
    
    % zonal accumulation
    fh(20) = figure(20),clf
    lat_plot=ybin(1:end-1)+.5;
    count_plot = nansum(COUNT,2);
    count_plot_wod_ctd = nansum(COUNT_CTD_WOD,2);
    count_plot_wod_pfl = nansum(COUNT_PFL_WOD,2);
    plot(lat_plot,count_plot,lat_plot,count_plot_wod_ctd,lat_plot,count_plot_wod_pfl)
    set(gca,'xlim',[-90 90],'ytick',0:5000:40000,'yticklabel',0:5000:40000)
    legend('MEOP-CTD','WOD13-CTD','WOD13-PFL')
    %title('Number of profiles per unit of latitude')
    xlabel('latitude'); ylabel('number of profiles (per unit of latitude)');
    nfile=sprintf('%sglobal/histogram',conf.mapsdir);
    print('-depsc2',nfile);
    print('-dpdf',nfile);
    print(nfile,'-dpng','-r300');
    %eval(['export_fig ' nfile '.png -m5 -nocrop']);
    %
    %
    % % total SH
    % figure(21);clf,
    % set(gcf, 'Color', 'w');
    % set(gcf, 'units','centimeters','Position', [2 5 20 20]);
    % m_proj('stereographic','lat',-90,'long',90,'radius',60);
    % m_grid('box','on','xtick',12,...
    %     'tickdir','out','ytick',[-90:10:-40],...
    %     'XaxisLocation','top','linest','-');
    % hold on
    % m_pcolor(XBIN,YBIN,log10(COUNT)); shading flat; caxis([0 3]);
    % topoplot([-180 180 -90 -30],0,2,'m_map')
    % colormap default;colorbar('ytick',0:3,'yticklabel',[1 10 100 1000])
    % suptitle(sprintf('MEOP-CTD SH dataset : %d profiles, %d deployments, %d tags',Nprof1,length(unique(Ndep1)),Ntag1));
    % nfile=sprintf('%sglobal/density_SH',conf.mapsdir);
    % eval(['export_fig ' nfile '.png -m5 -nocrop']);
    % eval(['export_fig ' nfile '.pdf -m5 -nocrop -painters']);
    %
    % % total NH
    % figure(22);clf,
    % set(gcf, 'Color', 'w');
    % set(gcf, 'units','centimeters','Position', [15 15 20 20]);
    % m_proj('stereographic','lat',90,'long',90,'radius',60);
    % m_grid('box','on','xtick',12,...
    %     'tickdir','out','ytick',[40:10:90],...
    %     'XaxisLocation','bottom','linest','-');
    % hold on
    % m_pcolor(XBIN,YBIN,log10(COUNT)); shading flat; caxis([0 3]);
    % topoplot([-180 180 30 90],0,2,'m_map')
    % colormap default;colorbar('ytick',0:3,'yticklabel',[1 10 100 1000])
    % suptitle(sprintf('MEOP-CTD NH dataset : %d profiles, %d deployments, %d tags',Nprof2,length(unique(Ndep2)),Ntag2));
    % nfile=sprintf('%sglobal/density_NH',conf.mapsdir);
    % eval(['export_fig ' nfile '.png -m5 -nocrop']);
    % eval(['export_fig ' nfile '.pdf -m5 -nocrop -painters']);
    %
    % total by nation
    fh(23) = figure(23);clf,
    set(gcf, 'Color', 'w');
    set(gcf, 'units','centimeters','Position', [5 30 20 14]);
    hold on
    pcolor(XBIN,YBIN,log10(COUNT)); shading flat; caxis([0 3]);
    topoplot([-250 110 -80 80],0,2); latratio(45)
    set(gca,'xtick',[-240:60:60],'xticklabel',{'120E','180E','120W','60W','0E','60E'});
    set(gca,'ytick',[-80:20:80],'yticklabel',{'80S','60S','40S','20S','0N','20N','40N','60N','80N'});
    grid on
    colormap default;colorbar('ytick',0:3,'yticklabel',[1 10 100 1000])
    suptitle(sprintf('MEOP-CTD dataset : %d profiles, %d deployments, %d tags',Nprof3,length(unique(Ndep3)),Ntag3));
    nfile=sprintf('%sglobal/density_global',conf.mapsdir);
    print(nfile,'-dpng','-r300');
    print(nfile,'-dpdf');

    %
    % % total WOD13-CTD
    % figure(24);clf,
    % set(gcf, 'Color', 'w');
    % set(gcf, 'units','centimeters','Position', [5 30 20 14]);
    % hold on
    % pcolor(XBIN,YBIN,log10(COUNT_CTD_WOD)); shading flat; caxis([0 3]);
    % topoplot([-250 110 -80 80],0,2); latratio(45)
    % set(gca,'xtick',[-240:60:60],'xticklabel',{'120E','180E','120W','60W','0E','60E'});
    % set(gca,'ytick',[-80:20:80],'yticklabel',{'80S','60S','40S','20S','0N','20N','40N','60N','80N'});
    % grid on
    % colormap default;colorbar('ytick',0:3,'yticklabel',[1 10 100 1000])
    % suptitle(sprintf('WOD13-CTD dataset : %d profiles',length(data_ctd_wod.latitude)));
    % nfile=sprintf('%sglobal/density_ctd_wod13',conf.mapsdir);
    % eval(['export_fig ' nfile '.png -m5 -nocrop']);
    % eval(['export_fig ' nfile '.pdf -m5 -nocrop -painters']);
    %
    % % total WOD13-PFL
    % figure(25),clf,
    % set(gcf, 'Color', 'w');
    % set(gcf, 'units','centimeters','Position', [5 30 20 14]);
    % hold on
    % pcolor(XBIN,YBIN,log10(COUNT_PFL_WOD)); shading flat; caxis([0 3]);
    % set(gca,'xtick',[-240:60:60],'xticklabel',{'120E','180E','120W','60W','0E','60E'});
    % set(gca,'ytick',[-80:20:80],'yticklabel',{'80S','60S','40S','20S','0N','20N','40N','60N','80N'});
    % grid on
    % topoplot([-250 110 -80 80],0,2); latratio(45)
    % colormap default;colorbar('ytick',0:3,'yticklabel',[1 10 100 1000])
    % suptitle(sprintf('WOD13-PFL dataset : %d profiles',length(data_pfl_wod.latitude)));
    % nfile=sprintf('%sglobal/density_pfl_wod13',conf.mapsdir);
    % eval(['export_fig ' nfile '.png -m5 -nocrop']);
    % eval(['export_fig ' nfile '.pdf -m5 -nocrop -painters']);
    
end

%%

% %% descriptive plots per nations
% for kk=1:length(conf.list_group),
%     group=conf.list_group{kk};
%     group=group(group~=' ');
%     Itag_tmp=conf.Itag_group{kk};
%     Mgroup=[];
%     t_grp=[];
%     for ii=1:length(Itag_tmp),
%         EXP = conf.lfic{Itag_tmp(ii)};
%         info_deployment=load_info_deployment(conf,EXP);
%         list_tag = dir([info_deployment.dir '*_prof.nc']);
%         if length(list_tag)>0
%             for jj=1:length(list_tag)
%                 name_prof = sprintf('%s%s',info_deployment.dir,list_tag(jj).name);
%                 M=ARGO_load_qc(name_prof,3);
%                 if length(M.JULD)>0
%                     Mgroup=ARGO_concat(Mgroup,M);
%                     GROUP_temp=cellstr(repmat(info_deployment.EXP,M.np,1));
%                     t_grp=[t_grp;GROUP_temp];
%                 end
% 
%             end
%             %N=Mgroup.np;
%             % Mgroup.GROUP=cellstr(repmat(info_deployment.EXP,N,1));
%         end
% 
%     end
%     if ~isempty(Mgroup)
%         Mgroup.GROUP=t_grp;
%         nfile=sprintf('%sdescriptive_%s.png',conf.maindir,group);
%         str_title=sprintf('%s data in MEOP-CTD',group);
%         figure(12),clf,descriptive_plot_ARGO(Mgroup,12,str_title,nfile,conf);
%     end
% end

% %% load all
% Mall=[];
% Mall_SO=[];
% Itag=[];
% Itag_SO=[];
% for kk=1:length(conf.list_group),
%     group=conf.list_group{kk};
%     Itag_tmp=conf.Itag_group{kk};
%     for ii=1:length(Itag_tmp),
%         Mgroup=[];        t_grp=[];
%         Mgroup_SO=[];     t_grp_SO=[];
%         EXP = conf.lfic{Itag_tmp(ii)};
%         info_deployment=load_info_deployment(conf,EXP);
%         list_tag = dir([info_deployment.dir '*_prof.nc']);
%         if length(list_tag)>0
%             for jj=1:length(list_tag)
%                 name_prof = sprintf('%s%s',info_deployment.dir,list_tag(jj).name);
%                 M=ARGO_load_qc(name_prof,3);
%                 if length(M.JULD)>0
%                     Mgroup=ARGO_concat(Mgroup,M);
%                     GROUP_temp=cellstr(repmat(info_deployment.NATION,M.np,1));
%                     t_grp=[t_grp;GROUP_temp];
%                     Itag=[Itag Itag_tmp(ii)];
%                     if nanmean(M.LATITUDE)<-20
%                         Mgroup_SO=ARGO_concat(Mgroup_SO,M);
%                         GROUP_temp=cellstr(repmat(info_deployment.NATION,M.np,1));
%                         t_grp_SO=[t_grp_SO;GROUP_temp];
%                         Itag_SO=[Itag_SO Itag_tmp(ii)];
%                     end
%                 end
%             end
%             if length(t_grp)>0
%                 Mgroup.GROUP=t_grp;
%                 Mall=ARGO_concat(Mall,Mgroup);
%             end
%             if length(t_grp_SO)>0
%                 Mgroup_SO.GROUP=t_grp_SO;
%                 Mall_SO=ARGO_concat(Mall_SO,Mgroup_SO);
%             end
%         end
% 
%     end
% end
% Itag=unique(Itag);
% Itag_SO=unique(Itag_SO);
% 
% %%
% 
% %
% nfile=sprintf('%sdescriptive.png',conf.maindir);
% str_title='MEOP-CTD dataset';
% figure(12),clf,descriptive_plot_ARGO(Mall,12,str_title,nfile);
% 
% %
% nfile=sprintf('%shistoARGO_SO.png',conf.maindir);
% figure(13),clf,histogram_ARGO(Mall_SO,13,nfile);
% feval(@print,'-dpng','-r300',[conf.maindir 'histoARGO_SO']);
% 



