function create_fr0(conf,EXP,one_smru_name)

if isempty(conf),
    conf = init_mirounga;
end

if ~exist('one_smru_name','var') % all tags from EXP deployment
    one_smru_name = '';
elseif isempty(EXP),
    EXP=EXP_from_smru_name(one_smru_name);
end

info_deployment = load_info_deployment(conf,EXP,one_smru_name);
PI=info_deployment.PI;
NATION=info_deployment.NATION;

% don't process it if no raw odv file
if ~exist([conf.rawdir info_deployment.nomfic]),
    fprintf('%s: no raw file. not processed.',EXP);
    return
end

diary_file = [info_deployment.dir EXP '_diary.txt'];
diary(diary_file)










% create netcdf HR traj and prof if FR data exist
list_tag = info_deployment.list_tag;
list_deployment_hr = conf.list_deployment_hr;
for index=1:length(list_tag)
    
    smru_name = info_deployment.list_smru_name{index};
    name_prof = sprintf('%s%s_lr0_prof.nc',info_deployment.dir,smru_name);
    if ~any(strcmp(list_deployment_hr.Properties.RowNames,smru_name)) || ~exist(name_prof,'file')
        continue
    end
    disp(['Process hr smru_name=' smru_name])
    
    % load reference locations
    [smru_prefix,Nsplit] = Nsplit_from_smru_name(smru_name);
    ptt = conf.platform_json(strcmp(conf.list_smru_platform_code,smru_prefix)).ptt;
    jul = ncread(name_prof,'JULD')+712224;
    locs = load_locs_data(conf,smru_name,ptt,jul);
    disp(sprintf('  SMRU %s, PTT %8s: use %s locations for fr data',smru_name,ptt,locs.type))
    
    % read parameter in list_deployment_hr
    num_file = list_deployment_hr{smru_name,'instr_id'};
    year     = list_deployment_hr{smru_name,'year'};
    prefix   = list_deployment_hr{smru_name,'prefix'};
    if ~isempty(prefix) & ~isnan(prefix),
        prefix = [num2str(prefix) '_'];
    else
        prefix = '';
    end
    continuous = list_deployment_hr{smru_name,'continuous'};
    
    name_hr_file     = sprintf('%s%d/%s%d_ctd.txt',conf.rawdir_hr,year,prefix,num_file);
    % Treat case where there is a numerical prefix: find it and add it to list_deployment_hr.csv
    if ~exist(name_hr_file,'file'),
        lfile = dir(sprintf('%s%d/*_%d_ctd.txt',conf.rawdir_hr,year,num_file));
        if ~isempty(lfile),
            name_hr_file = fullfile(lfile(1).folder,lfile(1).name);
            prefix = strrep(lfile(1).name,sprintf('_%d_ctd.txt',num_file),'');
            conf.list_deployment_hr{smru_name,'prefix'} = str2num(prefix);
            writetable(conf.list_deployment_hr,[conf.processdir 'list_deployment_hr.csv'],...
                'WriteRowNames',1,'Delimiter',',');
        else
            disp(sprintf('  %s not found',name_hr_file)); 
            continue,
        end
    end
    
    % load metadata in lr file
    data_att = ncloadatt_struct(name_prof);
    listatt = fieldnames(data_att);
    
    %% chargement des data du fichier ctd haute resolution
    hrdata = load_hr_data(name_hr_file,continuous);
    Nhr = length(hrdata.date);
    isfluo = hrdata.isfluo;
    isoxy = hrdata.isoxy;
    islight = hrdata.islight;
    
    % fonction de detection des plongees (identique a celle utilisee pour les tdr)
    tdr =[hrdata.date,hrdata.P,hrdata.T,hrdata.S,hrdata.F,hrdata.O,hrdata.L];
    if continuous,
        [statdives,info_ana_dives,statdivestxt,datadives,datadivestxt,chg,daindexes] = ...
            ana_dives_fabien(tdr);
        Ibeg = daindexes(2,:)';
        seuil_depth=2;% pour avoir les profils jusqu'à la surface
        Iend=[];
        for kk=1:length(statdives)-1
            K=find(tdr(chg(2,kk):chg(1,kk+1),2)<seuil_depth);
            if length(K)>0
               Iend(kk,1)=chg(2,kk)+K(1)-1; 
            else
                Iend(kk,1) = chg(2,kk);
            end
        
        end
        Iend(end+1)=chg(2,end);
    else
        Ibeg = [1;find(abs(diff(hrdata.P))>10)];
        Iend = [find(abs(diff(hrdata.P))>10)-1;length(Ibeg)];
        if Iend(1)==0,
            Ibeg(1)=[];
            Iend(1)=[];
        end
    end
    % remove last profile when it has a bad date
    I = find(diff(hrdata.date(Iend))<0);
    N = length(Ibeg);
    if length(I),
        Ibeg = Ibeg(1:I(1));
        Iend = Iend(1:I(1));
        disp(sprintf('  Date issue in %s hr data, %d instr_id: %d profiles removed',...
            smru_name,num_file,N-I(1)))
    end
    N = length(Ibeg);
    jul_hr = hrdata.date(Iend);       
        
    % associate locations
    try
        I = find(~isnan(locs.jul.*locs.lat.*locs.lon));
        lat_hr = interp1(locs.jul(I),locs.lat(I),jul_hr,'linear',NaN);
        lon_hr = interp1(locs.jul(I),locs.lon(I),jul_hr,'linear',NaN);
    catch
        disp(['  Problem locating hr profiles: ' name_hr_file])
    end
    
    if all(isnan(lat_hr))
        disp(['  Problem locating hr profiles: ' name_hr_file])
    end
    
    %% Save traj files (raw data)
    sc_write_traj_file;
    
    % creation des profils de remontees
    % PRES, TEMP, SALI, CHLA, DOXY, LIGHT, std_std --> **_std_lev
    PRES=NaN*zeros(6000,N);
    TEMP=NaN*zeros(6000,N);
    SALI=NaN*zeros(6000,N);
    CHLA=NaN*zeros(6000,N);
    DOXY=NaN*zeros(6000,N);
    LIGHT=NaN*zeros(6000,N);
    std_lev=[1:1000]';
    for ii=1:N
        dive=tdr(Ibeg(ii):Iend(ii),:);
        [m,J]=max(dive(:,2));
        dive(1:min(J+5,size(dive,1)),:)=[];
        PRES(1:size(dive,1),ii)=dive(:,2);
        TEMP(1:size(dive,1),ii)=dive(:,3);
        SALI(1:size(dive,1),ii)=dive(:,4);
        CHLA(1:size(dive,1),ii)=dive(:,5);
        DOXY(1:size(dive,1),ii)=dive(:,6);
        LIGHT(1:size(dive,1),ii)=dive(:,7);
    end
    sc_profiles_interp;
    P=PRES_std_lev;
    T=TEMP_std_lev;
    S=SALI_std_lev;
    F=CHLA_std_lev;
    O=DOXY_std_lev;
    L=LIGHT_std_lev;
     
    %% convert into fcell temporary format
    N = size(P,2);
    hi=zeros(N,10);
    hs={};
    for kk=1:N,
        hs(kk)={smru_name};
    end
    hi(:,2)=index;
    hi(:,3)=1:N;
    hi(:,4)=jul_hr;
    hi(:,5:6)=[lat_hr lon_hr];
    hi(:,9)=1;
    hi(:,10)=hi(:,2);
    
    PTi=cell(1,N); PSi=cell(1,N);
    PFi=cell(1,N); POi=cell(1,N); PLi=cell(1,N);
    T(T==999)=NaN; S(S==999)=NaN;
    F(F==999)=NaN; O(O==999)=NaN; L(L==999)=NaN;
    for ii=1:N,        
        PTi{ii}=[P(:,ii) T(:,ii)];
        PSi{ii}=[P(:,ii) S(:,ii)];
        PFi{ii}=[P(:,ii) F(:,ii)];
        POi{ii}=[P(:,ii) O(:,ii)];
        PLi{ii}=[P(:,ii) L(:,ii)];
        hi(ii,8)=sum(~isnan(T(:,ii)));        
    end
    
    % find bad location
    I=find(hi(:,5)~=0 & ~isnan(hi(:,4).*hi(:,5).*hi(:,6)));
    hi=hi(I,:);
    PTi=PTi(I); PSi=PSi(I);
    PFi=PFi(I); POi=POi(I); PLi=PLi(I);
    hs=hs(I);
    
    % correct dates: put them back to Matlab julian date norm
    hi(hi(:,4)<1e5,4)=hi(hi(:,4)<1e5,4)+693962;
    
    % sort date
    I=zeros(1,length(PTi));
    ltag=unique(hi(:,10));n=0;
    for kk=1:length(ltag),
        J=find(hi(:,10)==ltag(kk));
        [a,K]=sort(hi(J,4));
        I(n+1:n+length(J))=J(K);
        n=n+length(J);
    end
    hi=hi(I,:);
    PTi=PTi(I); PSi=PSi(I);
    PFi=PFi(I); POi=POi(I); PLi=PLi(I);
    hs=hs(I);
    
    %save fcell
    name_fcell=[conf.temporary_fcell smru_name '_fr0_fcell.mat'];
    save(name_fcell,'hi','hs','PTi','PSi',...
        'PFi','POi', 'PLi','EXP','PI','NATION','isoxy','isfluo','islight');
    
    %% save in Argo netcdf format
    if length(hi)>0
        suffix = 'fr0_prof.nc';
        convert_fcell2ARGO(conf,info_deployment.EXP,name_fcell,suffix,1000,smru_name);
    end
    disp(sprintf('  %s: hr tags with %d profiles',smru_name,N));
    
    % set default value for qc flags
    suffix = '_fr0';
    name_prof = sprintf('%s%s%s_prof.nc',info_deployment.dir,smru_name,suffix);
    sc_filtre_seals_qc;
    
end

diary off
