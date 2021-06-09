%% sc_plot_data_tags

if ~ismember(mode,{'raw','adj'}), 
    error('error: plotting mode must be one of raw/adj values'); 
end

close all
list_tag = getfield(info_deployment,['list_tag' suffix]);
plotdir = [conf.plotdir,EXP,'/'];
[s, mess, messid] = mkdir(plotdir);

for ii=1:length(list_tag),
    %%
    
    name_prof = sprintf('%s%s',info_deployment.dir,list_tag(ii).name);
    calib = ncload_struct(name_prof,'SCIENTIFIC_CALIB_COEFFICIENT');
    
    %descr=M.list_descr{ii}; TAG=descr(7:8);
    switch mode
        case 'raw'
            Mqc=ARGO_load_qc(name_prof,4);
        case 'adj'
            Mqc=ARGO_load_qc(name_prof,3);
        otherwise
            disp('error type mode');
    end
    if isfield (Mqc,'smru_platform_code')
        descr=Mqc.smru_platform_code';
    else
        descr=EXP';
    end
    Mqc.Tmask=double(Mqc.TEMP_QC<2); Mqc.TEMP(Mqc.TEMP_QC>2)=NaN;
    if isfield (Mqc,'PSAL')
        Mqc.Smask=double(Mqc.PSAL_QC<2); Mqc.PSAL(Mqc.PSAL_QC>2)=NaN;
        Ivalid = find(double(sum(Mqc.Tmask,1)~=0)+double(sum(Mqc.Smask,1)~=0));
    else
        Ivalid = find(double(sum(Mqc.Tmask,1)~=0));
        
    end
    
    if isempty(Ivalid), continue, end
    
    lat=Mqc.LATITUDE(Ivalid); lon=Mqc.LONGITUDE(Ivalid); is_lon_center_180=0;
    if any(lon<45&lon>-45),
        mlon=floor(min(lon)); Mlon=ceil(max(lon));
    else
        lon(lon<0)=lon(lon<0)+360;
        mlon=floor(min(lon)); Mlon=ceil(max(lon));
        is_lon_center_180=1;
    end
    lim=[floor(min(lat)) ceil(max(lat)) mlon Mlon];
    
    Mqc.PTMP = sw_ptmp ( Mqc.PSAL, Mqc.TEMP, Mqc.PRES, 0 );
    Mqc.SIG0 = sw_dens0( Mqc.PSAL, Mqc.PTMP ) - 1000 ;
    
    steps = [.01:.01:.05 .1:.1:.5 1:5];
    
    Tlim=[floor(min(Mqc.TEMP(:)-1)) ceil(max(Mqc.TEMP(:)+1))];
    [m,ii] = min(abs(steps-(Tlim(2)-Tlim(1))/10));
    Tcontour = Tlim(1):steps(ii):Tlim(2);
    if any(isnan(Tlim)), Tlim=[0 5]; Tcontour=0:5; end
    
    Slim=[floor(min(Mqc.PSAL(:)-0.1)) ceil(max(Mqc.PSAL(:)+0.1))];
    [m,ii] = min(abs(steps-(Slim(2)-Slim(1))/10));
    Scontour = Slim(1):steps(ii):Slim(2);
    if any(isnan(Slim)), Slim=[34 35]; Scontour=34:.1:35; end
    
    Dlim=[floor(min(Mqc.SIG0(:)-0.1)) ceil(max(Mqc.SIG0(:)+0.1))];
    [m,ii] = min(abs(steps-(Dlim(2)-Dlim(1))/10));
    Dcontour = Dlim(1):steps(ii):Dlim(2);
    if any(isnan(Dlim)), Dlim=[24 28]; Scontour=24:.5:28; end
    
    Fcontour=[0:8];
    Ocontour=[0:50:500];
    pmax = 1000;

    if isfield(plot_conf, 'lim'),  lim=plot_conf.lim; end
    if isfield(plot_conf,'pmax'), pmax=plot_conf.pmax; end
    if isfield(plot_conf,'Tlim'), Tlim=plot_conf.Tlim; end
    if isfield(plot_conf,'Slim'), Slim=plot_conf.Slim; end
    if isfield(plot_conf,'Dlim'), Dlim=plot_conf.Dlim; end
    if isfield(plot_conf,'Tcontour'), Tcontour=plot_conf.Tcontour; end
    if isfield(plot_conf,'Scontour'), Scontour=plot_conf.Scontour; end
    if isfield(plot_conf,'is_lon_center_180'), is_lon_center_180=plot_conf.is_lon_center_180; end
        
    % fiches
    Z=Mqc.JULD_LOCATION-min(Mqc.JULD_LOCATION);
    mZ=min(Z)-1; MZ=max(Z)+1; dZ=ceil((MZ-mZ)/10);
    edges=mZ:dZ:(MZ+dZ);
    %  smru_name=deblank(Mqc.smru_name{1});
    if isfield (Mqc,'smru_platform_code')
        smru_name=strrep(Mqc.smru_platform_code,'_','\_');
    else
        smru_name=strrep(Mqc.MEOP_CTD_NUMBER,'_','\_');
    end
    % calibration coefficient
    %I=find(strcmp(M.platform_number,M.list_descr{ii}));
    Tcoef = sscanf(calib.SCIENTIFIC_CALIB_COEFFICIENT(:,2,1,1)','t1= %f degC/km, t2= %f degC');
    if length(calib.SCIENTIFIC_CALIB_COEFFICIENT(1,:,1,1))>2
        Scoef = sscanf(calib.SCIENTIFIC_CALIB_COEFFICIENT(:,3,1,1)','s1= %f psu/km, s2= %f psu');
        Scoef2 = sscanf(calib.SCIENTIFIC_CALIB_COEFFICIENT(:,3,1,end)','s1= %f psu/km, s2= %f psu');
        isoffset = 0;
        if any(Scoef2-Scoef)
            for kk=2:size(calib.SCIENTIFIC_CALIB_COEFFICIENT,4),
                Scoef(:,kk) = sscanf(calib.SCIENTIFIC_CALIB_COEFFICIENT(:,3,1,kk)','s1= %f psu/km, s2= %f psu');
            end
            Scoef_mean = mean(Scoef,2);
            Scoef_std = std(Scoef,[],2);
            if any(isnan(Scoef_mean)), 
                Scoef_mean = [0 0]; 
                Scoef_std  = [0 0]; 
            end
            isoffset = 1;
        end
    end
    %%
    close all
    switch mode
        
        case 'raw'
            
            str={smru_name, ...
                'RAW DATA', ...
                sprintf('%d TEMP profiles',sum(double(sum(Mqc.Tmask,1)~=0)))};
            if isfield(Mqc,'Smask')
                str=[str sprintf('%d SALT profiles',sum(double(sum(Mqc.Smask,1)~=0)))];
            end
            str=[str sprintf('%d days',floor(max(Z)))];
            
        case 'adj'
            
            M=ARGO_load_qc(name_prof,2);
            str={smru_name, ...
                sprintf('%s data',mode), ...
                sprintf('%d TEMP profiles',sum(double(sum(Mqc.Tmask,1)~=0))), ...
                sprintf('%d days',floor(max(Z))),' ', ...
                sprintf('T1=%5.2f SI',Tcoef(1)),...
                sprintf('T0=%5.2f degC',Tcoef(2)),...
                sprintf('SigmaT= %5.2f',nanmean(M.TEMP_ADJUSTED_ERROR(:)))};
            
            if isfield(Mqc,'Smask')
                if isoffset
                    str=[str,sprintf('%d SALT profiles',sum(double(sum(Mqc.Smask,1)~=0))), ...
                        'variable salinity offset',...
                        sprintf('S1=%5.2f psu/km',Scoef_mean(1)),...
                        sprintf('S0=%5.2f +/- %5.2f psu',Scoef_mean(2),Scoef_std(2)),...
                        sprintf('SigmaS= %5.2f',nanmean(M.PSAL_ADJUSTED_ERROR(:)))];
                else
                    str=[str,sprintf('%d SALT profiles',sum(double(sum(Mqc.Smask,1)~=0))), ...
                        sprintf('S1=%5.2f psu/km',Scoef(1)),...
                        sprintf('S0=%5.2f psu',Scoef(2)),...
                        sprintf('SigmaS= %5.2f',nanmean(M.PSAL_ADJUSTED_ERROR(:)))];
                end
            end
            
        otherwise
            disp('error type mode');
            continue
            
    end
    
    nfile=sprintf('%s%s%s_diags_TS_%s.png',plotdir,descr,suffix,mode);
    if is_lon_center_180==1, Mqc.LONGITUDE(Mqc.LONGITUDE<0)=Mqc.LONGITUDE(Mqc.LONGITUDE<0)+360; end
    if ~isfield(Mqc,'PSAL')
        Mqc.PSAL=Mqc.TEMP*NaN;
    end
    TS_diags_ARGO(Mqc,Z,edges,str,lim,Tlim,Slim,Dlim,0,nfile);
    
    nfile=sprintf('%s%s%s_transect_%s.png',plotdir,descr,suffix,mode);
    transect_time_ARGO(Mqc,Tcontour,Scontour,edges,0,nfile,plot_conf);
    
    if isfield(Mqc,'CHLA'),
        nfile=sprintf('%s%s%s_transect_chla_%s.png',plotdir,descr,suffix,mode);
        transect_time_ARGO_CHLA(Mqc,Tcontour,Fcontour,edges,0,nfile,plot_conf);
    end
    
    if isfield(Mqc,'DOXY'),
        nfile=sprintf('%s%s%s_transect_doxy_%s.png',plotdir,descr,suffix,mode);
        transect_time_ARGO_DOXY(Mqc,Tcontour,Ocontour,edges,0,nfile,plot_conf);
    end
    
end
