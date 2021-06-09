function apply_tlc_fr(conf,EXP,one_smru_name)

if isempty(conf),
    conf = init_mirounga;
end

if ~exist('one_smru_name','var') % all tags from EXP deployment
    one_smru_name = '';
elseif isempty(EXP),
    EXP=EXP_from_smru_name(one_smru_name);
end


info_deployment=load_info_deployment(conf,EXP,one_smru_name);
if ~exist(info_deployment.dir), return, end

list_tag_fr0 = info_deployment.list_tag_fr0;
for ll=1:length(list_tag_fr0),
    
    % compute fr1 from fr0
    name_prof_fr0 = sprintf('%s%s',info_deployment.dir,info_deployment.list_tag_fr0(ll).name);
    Mqc=ARGO_load_qc(name_prof_fr0,0);
    suffix='_fr0'; smru_name = info_deployment.list_tag_fr0(ll).name(1:end-12);
    [Mqc,nn]=remove_profiles(info_deployment,...
        smru_name,'index',find(isnan(Mqc.LATITUDE.*Mqc.LONGITUDE.*Mqc.JULD)),suffix);
    ARGO_save_qc(name_prof_fr0,Mqc,0);
    
    name_prof_fr1 = strrep(name_prof_fr0,'_fr0_','_fr1_');
    if exist(name_prof_fr1,'file'), delete(name_prof_fr1); end
    copyfile(name_prof_fr0, name_prof_fr1,'f');
    Mqc = ARGO_load_qc(name_prof_fr1,1);
    
    % thermal cell correction
    [Tcorr1,Scorr1,Thp1] = thermal_cell_correction_2([.09 .05],.06,1,...
        Mqc.PSAL,Mqc.TEMP,Mqc.PRES,1);
    
    % density inversion
    c = gsw_SSO/35;
    CT = gsw_CT_from_t(Scorr1*c,Tcorr1,Mqc.PRES);
    Scorr2 = gsw_stabilise_SA_const_t(Scorr1*c,CT,Mqc.PRES)/c;
    
    % filter S
    X0 = Scorr2; [nr np] = size(X0);
    X1 = filter([1 4 6 4 1]/16,1,[nan(2,np); X0; nan(2,np)]);
    X2 = X1(5:end,:);
    X2(isnan(X2) & ~isnan(X0)) = X0(isnan(X2) & ~isnan(X0));
    Mqc.PSAL_ADJUSTED = X2;
    
    % filter T
    X0 = Tcorr1; [nr np] = size(X0);
    X1 = filter([1 4 6 4 1]/16,1,[nan(2,np); X0; nan(2,np)]);
    X2 = X1(5:end,:);
    X2(isnan(X2) & ~isnan(X0)) = X0(isnan(X2) & ~isnan(X0));
    Mqc.TEMP_ADJUSTED = X2;
    
    ncwrite(name_prof_fr1,'TEMP_ADJUSTED',Mqc.TEMP_ADJUSTED);
    ncwrite(name_prof_fr1,'PSAL_ADJUSTED',Mqc.PSAL_ADJUSTED);
    ncwriteatt(name_prof_fr1,'/','thermal_lag_adjustment','yes');
    
end

