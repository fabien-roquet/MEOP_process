function Mgroup=load_deployment(list)
list_tag=[];
ind=1;
if ~iscell(list)
    if (isnumeric(list) && size(list,2)==1)
        EXP = conf.lfic{list};
        info_deployment=load_info_deployment(conf,EXP);
        for ii=1:length(info_deployment.list_tag)
            list_tag{ind,1}=[info_deployment.dir,info_deployment.list_tag(jj).name];
            ind=ind+1;
        end
    elseif isnumeric(list)
        for ii=1:length(list)
            EXP = conf.lfic{list(ii)};
            info_deployment=load_info_deployment(conf,EXP);
            for jj=1:length(info_deployment.list_tag)
                list_tag{ind,1}=[info_deployment.dir,info_deployment.list_tag(jj).name];
                ind=ind+1;
            end
        end
    elseif ischar(list)
        list=find(strcmp(conf.lfic,list));
        EXP = conf.lfic{list};
        info_deployment=load_info_deployment(conf,EXP);
        for ii=1:length(info_deployment.list_tag)
            list_tag{ind,1}=[info_deployment.dir,info_deployment.list_tag(jj).name];
            ind=ind+1;
        end
    end
    
elseif iscell(list)
    for ii=1:length(list)
        tmp=find(strcmp(conf.lfic,list(ii)));
        EXP = conf.lfic{tmp(1)};
        info_deployment=load_info_deployment(conf,EXP);
        for jj=1:length(info_deployment.list_tag)
            list_tag{ind,1}=[info_deployment.dir,info_deployment.list_tag(jj).name];
            ind=ind+1;
        end
    end
end

Mgroup=[]; list_smru_platform={};
for itag=1:length(list_tag)
    name_prof = list_tag{itag};
    if exist(name_prof,'file'),
        Mqc=ARGO_load_qc(name_prof,1);
        list_smru_platform{itag}=Mqc.smru_platform_code;
        Mgroup=ARGO_concat(Mgroup,Mqc);
    else
        disp(['no data file: ' name_prof]);
        info_deployment.list_tag(itag)=[];
        calib_coeff(itag,:)=[];
    end
end
Mgroup.list_smru_platform=list_smru_platform;