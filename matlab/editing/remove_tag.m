function argo_qc = remove_tag (info_deployment,smru_name)
% remove a list of tags

smru_name = strrep(smru_name,'_lr0','');

name_prof = sprintf('%s%s%s_prof.nc',info_deployment.dir,smru_name,'_lr0');

if ~exist(name_prof,'file'), return; end

argo_qc=ARGO_load_qc(name_prof,0);
if any(argo_qc.PRES_QC(:)<=1)
    argo_qc.PRES_QC(argo_qc.PRES_QC~=9)=4;
    argo_qc.TEMP_QC(argo_qc.TEMP_QC~=9)=4;
    argo_qc.PSAL_QC(argo_qc.PSAL_QC~=9)=4;
    if isfield(argo_qc,'CHLA'),
        argo_qc.CHLA_QC(argo_qc.CHLA_QC~=9)=4;
    end
    if isfield(argo_qc,'DOXY'),
        argo_qc.DOXY_QC(argo_qc.DOXY_QC~=9)=4;
    end
    disp(sprintf('Tag %s lr0 removed',smru_name));
    ARGO_save_qc(name_prof,argo_qc,0);    
end

name_prof = sprintf('%s%s%s_prof.nc',info_deployment.dir,smru_name,'_fr0');
if ~exist(name_prof,'file'), return; end
argo_qc=ARGO_load_qc(name_prof,0);
if any(argo_qc.PRES_QC(:)<=1)
    argo_qc.PRES_QC(argo_qc.PRES_QC~=9)=4;
    argo_qc.TEMP_QC(argo_qc.TEMP_QC~=9)=4;
    argo_qc.PSAL_QC(argo_qc.PSAL_QC~=9)=4;
    if isfield(argo_qc,'CHLA'),
        argo_qc.CHLA_QC(argo_qc.CHLA_QC~=9)=4;
    end
    if isfield(argo_qc,'DOXY'),
        argo_qc.DOXY_QC(argo_qc.DOXY_QC~=9)=4;
    end
    
    disp(sprintf('Tag %s fr0 removed',smru_name));
    ARGO_save_qc(name_prof,argo_qc,0);
end
