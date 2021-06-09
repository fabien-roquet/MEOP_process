function locs = load_locs_data(conf,smru_name,ptt,jul)

    if isempty(conf),
        conf = init_mirounga;
    end

    [smru_prefix,Nsplit] = Nsplit_from_smru_name(smru_name);
    locs = load_crawl_data(conf,smru_prefix,ptt,jul);
    if isempty(locs), 
        locs = load_cls_data(conf,ptt,jul);
        if isempty(locs),
            info_deployment = load_info_deployment(conf,EXP_from_smru_name(smru_name),smru_name);
            name_prof = sprintf('%s%s_lr0_prof.nc',info_deployment.dir,smru_name);
            if exist(name_prof,'file')
                locs.type='smru';
                Mqc = ARGO_load_qc(name_prof,1);
                locs.jul = Mqc.JULD;
                locs.lat = Mqc.LATITUDE;
                locs.lon = Mqc.LONGITUDE;
                % remove doublons and sort by chronological order
                [jul,ij]=unique(locs.jul);
                locs.jul=jul;
                locs.lat=locs.lat(ij);
                locs.lon=locs.lon(ij);
            else
                locs = [];
                locs.type='none';
                locs.jul=[];
                locs.lat=[];
                locs.lon=[];
            end
        else
            locs.type='cls';
        end
    else
        locs.type='crawl';        
    end
    


function locs = load_cls_data(conf,ptt,jul)

    locs=[];

    datemin = min(jul);
    datemax = max(jul);
    I = find(ismember(conf.cls.ptt,ptt));
    if isempty(I),
        return
    end

    J=[];
    for ii=1:length(I),
        if datemin > conf.cls.datemin_jul(I(ii))-300 ...
                && datemax < conf.cls.datemax_jul(I(ii))+300
            J=I(ii);
        end
    end
    if isempty(J), 
        return
    end
    I=J;

    cls_file_name = conf.cls.list(I).name;
    fid = fopen([conf.cls.locdir,cls_file_name]);
    data = textscan(fid,...
        '%*s%f%f%f%f%s%s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s%*s',...
        'delimiter',';','headerlines',1);
    fclose(fid);
    clslat = data{1};
    clslon = data{2}; clslon(clslon>180)=clslon(clslon>180)-360;
    clslat2= data{3};
    clslon2= data{4}; clslon2(clslon2>180)=clslon2(clslon2>180)-360;
    clsqua = data{5};
    clsjul = datenum(data{6});
    [a,J]  = sort(clsjul);
    [c,ia,ic] = unique(clsjul(J)); J=J(ia);
    ia = find(~strcmp(clsqua(J),'Z')); J=J(ia);
    diff1 = abs(clslon2-clslon)+abs(clslat2-clslat);
    I=find(diff1~=0); I(I==1)=[];
    dist1= abs(clslon(I)-clslon(I-1))+abs(clslat(I)-clslat(I-1));
    dist2= abs(clslon2(I)-clslon2(I-1))+abs(clslat2(I)-clslat2(I-1));
    for ii=1:length(I),
        if dist2(ii)<dist1(ii)
            clslon(I(ii))=clslon2(I(ii));
            clslat(I(ii))=clslat2(I(ii));
        end
    end
    %figure(2),clf,plot(clslon,clslat,'+-',clslon2,clslat2,'+-')

    I=find(diff(clsjul)>0);
    locs.jul = clsjul(I);
    locs.lat = clslat(I);
    locs.lon = clslon(I);


function locs = load_crawl_data(conf,smru_name,ptt,jul)
    warning('off','MATLAB:table:ModifiedAndSavedVarnames')
    locs=[];

    datemin = min(jul);
    datemax = max(jul);
    I = find(ismember(conf.crawl.ptt,ptt) & ismember(conf.crawl.smru_name,smru_name));
    if isempty(I), 
        return
    end

    file_name = fullfile(conf.crawl.locdir,conf.crawl.list(I(1)).name);
    opts = detectImportOptions(file_name,'delimiter',',','ReadVariableNames',1);
    opts.SelectedVariableNames = {'GMT','mu_x','mu_y'};
    data = readtable(file_name,opts);
    data.Properties.VariableNames = {'date','lon','lat'};
    data{data{:,'lon'}>180,'lon'}=data{data{:,'lon'}>180,'lon'}-360;

    locs.jul = datenum(data{:,'date'});
    locs.lat = data{:,'lat'};
    locs.lon = data{:,'lon'};

    I=find(locs.jul>datemin-10 & locs.jul<datemax+10);
    if isempty(I)
        locs=[];
    else
        locs.jul = locs.jul(I);
        locs.lat = locs.lat(I);
        locs.lon = locs.lon(I);
    end
    

