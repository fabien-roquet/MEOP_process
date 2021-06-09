%% sc_load_odv_fcell_argo

EXP=info_deployment.EXP;
PI=info_deployment.PI;
NATION=info_deployment.NATION;

name_fcell1=[conf.rawdir info_deployment.EXP '_fcell.mat'];
load(name_fcell1);
name_fcell2=[conf.temporary_fcell info_deployment.EXP '_lr0_fcell.mat'];
copyfile(name_fcell1,name_fcell2);
save(name_fcell2,'NATION','-append');

N=length(PTi);
ltag=unique(hs);

%% save in Argo netcdf format
if length(hi)>0
    suffix = 'lr0_prof.nc';
    convert_fcell2ARGO(conf,info_deployment.EXP,name_fcell2,suffix,[],one_smru_name);
end

disp(sprintf('\t%d tags',length(ltag)));
disp(sprintf('\t%d profiles',N));


