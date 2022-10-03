function apply_tlc(conf,EXP,one_smru_name)

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

list_tag_hr0 = info_deployment.list_tag_hr0;
for ll=1:length(list_tag_hr0),
    
    % compute hr1 from hr0
    name_prof_hr0 = sprintf('%s%s',info_deployment.dir,info_deployment.list_tag_hr0(ll).name);
    name_prof_hr1 = strrep(name_prof_hr0,'_hr0_','_hr1_');
    if exist(name_prof_hr1,'file'), delete(name_prof_hr1); end
    copyfile(name_prof_hr0, name_prof_hr1);
    Mqc = ARGO_load_qc(name_prof_hr0,1);
    
    % thermal cell correction
    %[Tcorr1,Scorr1,Thp1] = thermal_cell_correction_2([.09 .05],.06,1,...
    %    Mqc.PSAL,Mqc.TEMP,Mqc.PRES,1);
    
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
    
    ncwrite(name_prof_hr1,'TEMP_ADJUSTED',Mqc.TEMP_ADJUSTED);
    ncwrite(name_prof_hr1,'PSAL_ADJUSTED',Mqc.PSAL_ADJUSTED);
    ncwriteatt(name_prof_hr1,'/','thermal_lag_adjustment','no');
    
    % compute lr1 from hr1
    name_prof_hr0 = sprintf('%s%s',info_deployment.dir,info_deployment.list_tag_hr0(ll).name);
    name_prof_hr1 = strrep(name_prof_hr0,'_hr0_','_hr1_');
    name_prof_lr0 = strrep(name_prof_hr0,'_hr0_','_lr0_');
    name_prof_lr1 = strrep(name_prof_hr0,'_hr0_','_lr1_');
    if exist(name_prof_lr1,'file'), delete(name_prof_lr1); end
    copyfile(name_prof_lr0, name_prof_lr1);
    Mqc_hr = ARGO_load_qc(name_prof_hr1,1);
    Mqc_lr = ARGO_load_qc(name_prof_lr1,1);
    Mqc_lr0 = ARGO_load_qc(name_prof_lr0,1);
    
    P_lr = Mqc_lr.PRES;
    T_lr = P_lr*NaN;
    S_lr = P_lr*NaN;
    [nr np] = size(P_lr);
    
    P_hr = Mqc_hr.PRES;
    T_hr = Mqc_hr.TEMP;
    S_hr = Mqc_hr.PSAL;
    
    for kk=1:np,
        
        I = find(~isnan(P_hr(:,kk).*T_hr(:,kk))); Ti = P_lr(:,kk)*NaN;
        if length(I)>1, Ti = interp1(P_hr(I,kk),T_hr(I,kk),P_lr(:,kk),'nearest',NaN); end
        T_lr(:,kk) = Ti;
        
        I = find(~isnan(P_hr(:,kk).*S_hr(:,kk))); Si = P_lr(:,kk)*NaN;
        if length(I)>1, Si = interp1(P_hr(I,kk),S_hr(I,kk),P_lr(:,kk),'nearest',NaN); end
        S_lr(:,kk) = Si;
        
    end
    
    T_lr(isnan(T_lr)) = Mqc_lr.TEMP(isnan(T_lr));
    S_lr(isnan(S_lr)) = Mqc_lr.PSAL(isnan(S_lr));
    
    ncwrite(name_prof_lr1,'TEMP_ADJUSTED',T_lr);
    ncwrite(name_prof_lr1,'PSAL_ADJUSTED',S_lr);
    ncwriteatt(name_prof_lr1,'/','thermal_lag_adjustment','no');
    
end


