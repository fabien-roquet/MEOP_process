%% script for adjustment of locations based on CLS post-processed data

%%
info_deployment = load_info_deployment(conf,EXP,one_smru_name);
list_tag = info_deployment.list_tag;
for index=1:length(list_tag),
    
    smru_name = info_deployment.list_smru_name{index};
    name_prof = sprintf('%s%s%s_prof.nc',info_deployment.dir,smru_name,'_lr0');
    ptt = ncreadatt (name_prof,'/','ptt');
    jul = ncread(name_prof,'JULD')+712224;
    locsystem = ncreadatt (name_prof,'/','loc_algorithm');
    plot_on=1; if strcmp(locsystem,'S'), plot_on=0; end
    locs = load_locs_data(conf,smru_name,ptt,jul);
    if strcmp(locs.type,'none'),
        continue
    end
    
    suffix = '_lr0';
    name_prof = sprintf('%s%s%s_prof.nc',info_deployment.dir,smru_name,suffix);
    if ~exist(name_prof,'file'), continue; end
    
    finfo = ncinfo(name_prof);
    dimNames = {finfo.Dimensions.Name};
    dimMatch = strcmp(dimNames,'N_PROF');
    Nprof = finfo.Dimensions(dimMatch).Length;
    lat = ncread(name_prof,'LATITUDE');
    lon = ncread(name_prof,'LONGITUDE');
    jul = ncread(name_prof,'JULD')+712224;
    ktag = find(ismember(conf.list_deployment_code,EXP));
    if strcmp(conf.platform_json(ktag).loc_algorithm,'K')
        ncwriteatt(name_prof,'/','loc_algorithm','CLS KALMAN');
        ncwrite   (name_prof,'POSITIONING_SYSTEM',repmat('K       ',Nprof,1)');
    end
    if strcmp(conf.platform_json(ktag).loc_algorithm,'L')
        ncwriteatt(name_prof,'/','loc_algorithm','CLS LEAST SQUARES');
        ncwrite   (name_prof,'POSITIONING_SYSTEM',repmat('LS      ',Nprof,1)');
    end
    if length(lat)<3, continue, end 

    newlat = interp1(locs.jul,locs.lat,jul,'linear','extrap');
    newlon = interp1(locs.jul,locs.lon,jul,'linear','extrap');
    Inan = find(jul < min(locs.jul)-2 | jul > max(locs.jul)+2);
    newlat(Inan) = lat(Inan);
    newlon(Inan) = lon(Inan);

    % speed filter 3m/s
    dist = spheredist(newlat,newlon);
    velo = diff(dist)./diff(jul)/86.4; % in m/s
    I=find(velo>3); Ibad = intersect(I,I+1); Ibad=setdiff(Ibad,[1 length(dist)]);
    if length(Ibad)>0
        jul1=jul; jul1(Ibad)=[]; 
        lat1=newlat; lat1(Ibad)=[]; 
        lon1=newlon; lon1(Ibad)=[];
        I = find(diff(jul1)==0);
        if length(I)
            jul1(I)=[]; lat1(I)=[]; lon1(I)=[];
        end
        newlat(Ibad) = interp1(jul1,lat1,jul(Ibad));
        newlon(Ibad) = interp1(jul1,lon1,jul(Ibad));
    end

    switch locs.type
        case 'cls'
            ncwriteatt(name_prof,'/','loc_algorithm','CLS SMOOTH KALMAN');
            disp(sprintf('SMRU %s, PTT %8s: use %s locations',smru_name,ptt,locs.type))
        case 'crawl'
            ncwriteatt(name_prof,'/','loc_algorithm','CRAWL');
            disp(sprintf('SMRU %s, PTT %8s: use %s locations',smru_name,ptt,locs.type))
    end
    ncwrite   (name_prof,'LATITUDE' ,newlat);
    ncwrite   (name_prof,'LONGITUDE',newlon);

    if plot_on
        figure('visible','off') ,clf
        plot(locs.lon,locs.lat,'-'),hold on
        xlim(get(gca,'xlim'))
        ylim(get(gca,'ylim'))
        plot(lon,lat,'x',newlon,newlat,'+')
        title(sprintf('%s loc: %8s %s',locs.type,ptt,smru_name))
        legend(locs.type,'old lr','new lr')
        format_figure_centred([15 15])
        plotdir = [conf.plotdir,info_deployment.EXP,'/'];
        [s,mess,messid] = mkdir(plotdir);
        print(sprintf('%slocation_%s_%08s',plotdir,smru_name,ptt),'-dpng','-r300')
    end
    
end

