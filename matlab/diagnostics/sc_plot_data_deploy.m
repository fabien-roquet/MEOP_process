%% sc_plot_data_deploy

close all;
list_tag = dir([info_deployment.dir '*' suffix '_prof.nc']);
plotdir = [conf.plotdir,EXP,'/'];
[s, mess, messid] = mkdir(plotdir);

load_qc=0;
if strcmp(mode,'adj'), load_qc=3; end

Mqc=[];
for ii=1:length(list_tag)
    name_prof = sprintf('%s%s',info_deployment.dir,list_tag(ii).name);
    arg=ARGO_load_qc(name_prof,load_qc);
    if ~isfield(arg,'PSAL')
        arg.PSAL=arg.TEMP*NaN;
        arg.PSAL_QC=arg.TEMP_QC*0+9;
    end
    arg.index_tag(:)= ii;
    if length(arg.JULD)>0
        Mqc=ARGO_concat(Mqc,arg);
    end
end

if ~isempty(Mqc)
    
    if length(Mqc.index_tag)>length(Mqc.JULD)
        Mqc.index_tag(length(Mqc.JULD)+1:end)=[];
    end
    
    Mqc.Tmask=double(Mqc.TEMP_QC<2); Mqc.TEMP(Mqc.TEMP_QC>2)=NaN;
    Mqc.Smask=double(Mqc.PSAL_QC<2); Mqc.PSAL(Mqc.PSAL_QC>2)=NaN;
    
    Ivalid = find(double(sum(Mqc.Tmask,1)~=0)+double(sum(Mqc.Smask,1)~=0));
    
    if ~isempty(Ivalid),
        
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
        
        Z=Mqc.index_tag; Ncol=max(Z); mZ=min(Z); MZ=max(Z); edges=mZ:MZ;
        str= { EXP,' ', ...
            sprintf('%d SRDL tags',length(Mqc.list_descr)), ...
            sprintf('%d TEMP profiles',sum(double(sum(Mqc.Tmask)~=0))), ...
            sprintf('%d SALT profiles',sum(double(sum(Mqc.Smask)~=0))), ...
            sprintf('%d T+S  profiles',sum(double(sum(Mqc.Tmask.*Mqc.Smask)~=0))) };
        if is_lon_center_180==1,
            Mqc.LONGITUDE(Mqc.LONGITUDE<0)=Mqc.LONGITUDE(Mqc.LONGITUDE<0)+360;
        end
        
        printname=[plotdir EXP suffix '_recapARGO_' mode '.png'];
        nfile=sprintf('%s%s%s_histoARGO_%s.png',plotdir,EXP,suffix,mode);
        TS_diags_ARGO(Mqc,Z,edges,str,lim,Tlim,Slim,Dlim,0,printname);
        histogram_ARGO(Mqc,0,nfile);
        
    end
    
end

