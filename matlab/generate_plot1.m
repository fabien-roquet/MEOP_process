function generate_plot1(conf,EXP,one_smru_name)


if isstr(conf),
    plot1_mode = conf;
    conf = init_mirounga;
    conf.plot1_mode = plot1_mode;
elseif isempty(conf),
    conf = init_mirounga;
end

if ~exist('one_smru_name','var') % all tags from EXP deployment
    one_smru_name = '';
elseif isempty(EXP),
    EXP=EXP_from_smru_name(one_smru_name);
end

if ~isfield(conf,'plot1_mode')
    conf.plot1_mode = 'fast';
end

close all

if ~exist('one_smru_name','var') % all tags from EXP deployment
    info_deployment=load_info_deployment(conf,EXP);
    if isempty(info_deployment.list_smru_name)
        return
    end
    [s,mess,messid] = mkdir(sprintf('%s%s',conf.calibplotdir,info_deployment.EXP));
    delete(sprintf('%s%s/*.png',conf.calibplotdir,info_deployment.EXP))
    disp(['calibration plots: ' EXP]);
else  % tag smru_tag only
    info_deployment=load_info_deployment(conf,EXP,one_smru_name);
    if isempty(info_deployment.list_smru_name)
        return
    end
    [s,mess,messid] = mkdir(sprintf('%s%s',conf.calibplotdir,info_deployment.EXP));
    delete(sprintf('%s%s/calibration_%s_*.png',conf.calibplotdir,info_deployment.EXP,one_smru_name))
    disp(['calibration plots: ' one_smru_name]);
end 


%% visualization of TS plots
conf_clim.lon=[];
conf_clim.lat=[];
for ii = 1:length(info_deployment.list_tag)
    conf_clim.lon=[conf_clim.lon;ncread([info_deployment.dir info_deployment.list_tag(ii).name],'LONGITUDE')];
    conf_clim.lat=[conf_clim.lat;ncread([info_deployment.dir info_deployment.list_tag(ii).name],'LATITUDE')];
end

conf_clim.woddir=conf.woddir;
conf_clim.coradir=conf.coradir;
conf_clim.meopdir=conf.meopdir;

argo_wod  = []; %load_WOD(conf_clim);
argo_cora = load_CORA(conf_clim);
argo_meop = load_MEOP(conf_clim);


%% load seal deployment data
Mgroup=[]; list_smru_platform={};list_platform_number={};
for itag = 1:length(info_deployment.list_tag)
    name_prof = sprintf('%s%s',info_deployment.dir,...
        strrep(info_deployment.list_tag(itag).name,'_lr0','_lr1'));
    if exist(name_prof,'file'),
        Mqc=ARGO_load_qc(name_prof,3);
        list_smru_platform{itag}=Mqc.smru_platform_code;
        if length(Mqc.platform_number),
            list_platform_number{itag}=Mqc.platform_number{1};
        else
            list_platform_number{itag}='';
        end
        Mgroup=ARGO_concat(Mgroup,Mqc);
    else
        disp(['no data file: ' name_prof]);
        info_deployment.list_tag(itag)=[];
    end
end
Mgroup.list_smru_platform=list_smru_platform;
Mgroup.list_platform_number=list_platform_number;

%%
[I,J]=find(Mgroup.TEMP_QC==1);
J=unique(J);
if length(J)>0
    tmp=Mgroup.TEMP(:,J);
    tmp=nanmax(tmp);
    tmp=sort(tmp(~isnan(tmp)));
    if length(tmp)>0
        temp_max=tmp(min([ceil(.99*length(tmp)) length(tmp)]));
        tmp=Mgroup.TEMP(:,J);
        tmp=nanmin(tmp);
        tmp=sort(tmp(~isnan(tmp)));
        temp_min=tmp(max([floor(.01*length(tmp)) 1]));
    else
        temp_max=5;
        temp_min=-2;
    end
else
    temp_max=5;
    temp_min=-2;
end
[I,J]=find(Mgroup.PSAL_QC==1);
J=unique(J);
if length(J)>0
    tmp=Mgroup.PSAL(:,J);
    tmp=nanmax(tmp);
    tmp=sort(tmp(~isnan(tmp)));
    if length(tmp)>0
        sal_max=max([35,tmp(min([ceil(.99*length(tmp)) length(tmp)]))]);
        tmp=Mgroup.PSAL(:,J);
        tmp=nanmin(tmp);
        tmp=sort(tmp(~isnan(tmp)));
        sal_min=tmp(max([floor(.01*length(tmp)) 1]));
    else
        sal_max=36;
        sal_min=33;
    end
else
    sal_max=36;
    sal_min=33;
end

%% compute comparison plots
for itag = 1:length(info_deployment.list_tag)
    
    %%
    conf_adjustement=[];
    conf_adjustement.Nprof_diags = 200;
    
    I=find(strcmp(Mgroup.platform_number,Mgroup.list_platform_number{itag}));
    if length(I)==length(Mgroup.index_tag)
        Mtemp=Mgroup;
    elseif ~isempty(I)
        Mtemp=extract_profil(Mgroup,I);
    else
        continue
    end
    conf_adjustement.N_profiles = length(I);
    conf_adjustement.argo_qc = Mtemp;
    
    coeff = conf.table_coeff;
    salinity_offsets = conf.table_salinity_offsets;
    
    smru_name = info_deployment.list_smru_name{itag};
    P1=0; P2=0; T1=0; T2=0; S1=0; S2=0; F1=0.6; F2=0;
    list_var = {'T1','T2','S1','S2'};
    for kk = 1:length(list_var),
        if any(strcmp(list_var{kk},coeff(smru_name,:).Properties.VariableNames)) & ...
                coeff{smru_name,list_var{kk}} & ~isnan(coeff{smru_name,list_var{kk}})
            eval([list_var{kk} ' = coeff{smru_name,list_var{kk}};'])
        end
    end
    conf_adjustement.P1 = P1;
    conf_adjustement.P2 = P2;
    conf_adjustement.T1 = T1;
    conf_adjustement.T2 = T2;
    conf_adjustement.S1 = S1;
    conf_adjustement.S2 = S2;
    conf_adjustement.F1 = F1;
    conf_adjustement.F2 = F2;
    conf_adjustement.offset = load_salinity_offset(smru_name,salinity_offsets,conf_adjustement.N_profiles);
    
    conf_adjustement.pause=0;
    conf_adjustement.Tlim=[temp_min temp_max]; 
    conf_adjustement.Slim=[sal_min sal_max];
    conf_adjustement.hfig = 0;
    
    conf_adjustement.nomfig=sprintf('%s/%s/calibration_%s_0.png',...
        conf.calibplotdir,info_deployment.EXP,smru_name);
    
    if ~isempty(argo_wod), conf_adjustement.argo_wod = argo_wod; end
    if ~isempty(argo_cora), conf_adjustement.argo_cora = argo_cora; end
    if ~isempty(argo_meop), conf_adjustement.argo_meop= argo_meop; end
    
    conf_adjustement.argo_qc2 = Mgroup;
    conf_adjustement.nomfig2=sprintf('%s/%s/calibration_%s_other_tags.png',...
        conf.calibplotdir,info_deployment.EXP,smru_name);
    TS_diags_comparison(conf_adjustement);
    
    
    if ~strcmp(conf.plot1_mode,'fast')
        conf_adjustement=rmfield(conf_adjustement,'argo_qc2');
        J=1:conf_adjustement.Nprof_diags:conf_adjustement.N_profiles;
        if J(end)~=conf_adjustement.N_profiles
            J(end+1)=conf_adjustement.N_profiles;
        end
        for jj=1:length(J)-1
            Msmall = extract_profil(Mtemp,J(jj):J(jj+1)-1);
            if conf_adjustement.offset
                conf_adjustement.offset=Offset(J(jj):J(jj+1)-1);
            end
            conf_adjustement.argo_qc = Msmall;
            conf_adjustement.nomfig=sprintf('%s%s/calibration_%s_%03d_part%02d.png',...
                conf.calibplotdir,info_deployment.EXP,...
                smru_name,conf_adjustement.Nprof_diags,jj);
            TS_diags_comparison(conf_adjustement);
        end
    end
    
    %%
    % %% estimate Toffset using freezing temperature
    % Mqc=ARGO_load_qc(name_prof,1);
    % Toffset = estimate_Toffset_freezing_temp(Mqc);
    %
    % %% estimate offset minimizing surface salinity differences
    % %  (additive offset, put in salinity_offset)
    %
    % Mqc=ARGO_load_qc(name_prof,1);
    % Soffset = estimate_Soffset_surface_salinity(Mqc,'NKP');
    %
    %
    
end



