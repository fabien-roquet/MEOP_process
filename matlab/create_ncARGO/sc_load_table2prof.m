%% sc_load_odv_fcell_argo

EXP=info_deployment.EXP;
PI=info_deployment.PI;
NATION=info_deployment.NATION;

% read header
clear F O C L;
isfluo=0; isoxy=0; iscond=0;islight=0;
%Fluo=strfind(tline,'Fluorescence');
%Oxy =strfind(tline,'Oxygen');
%Cond=strfind(tline,'Conductivity');
%Light=strfind(tline,'Light');

% read data
smru_names = table.smru_name;
N = length(table.smru_name);
numProf = (1:N)';
Ihead=find(numProf~=0);
date = table.END_DATE;
lon = table.LON;
lat = table.LAT;

% create header hi (fcell format)
hi=zeros(N,10);
hs=smru_names(Ihead);
ltag=unique(hs);
for ii=1:length(ltag),
    hi(strcmp(ltag{ii},hs),2)=ii;
end
hi(:,3)=numProf(Ihead); % not used
hi(:,4)=datenum(date(Ihead),'yyyy-mm-dd HH:MM');
hi(:,5:6)=[lat(Ihead) lon(Ihead)];
hi(:,9)=1;
hi(:,10)=hi(:,2);



PTi=cell(1,N); PSi=cell(1,N); 
for kprof = 1: N,

    str = table.TEMP_DBAR{kprof}; str = strrep(str,',',' ');
    Pt = sscanf(str,'%f');

    str = table.TEMP_VALS{kprof}; str = strrep(str,',',' ');
    Tt = sscanf(str,'%f');

    str = table.SAL_DBAR{kprof}; str = strrep(str,',',' ');
    Ps = sscanf(str,'%f');

    str = table.SAL_VALS{kprof}; str = strrep(str,',',' ');
    Ss = sscanf(str,'%f');

    if table.N_TEMP(kprof)>0,
        PTi{kprof}=[Pt Tt];
    else
        PTi{kprof}=zeros(8,2)*NaN;
    end
    
    if table.N_SAL(kprof)>0,
        PSi{kprof}=[Ps Ss];
    else
        PSi{kprof}=zeros(8,2)*NaN;
    end
    
    PFi{kprof}=PTi{kprof}*NaN;
    POi{kprof}=PTi{kprof}*NaN;
    PLi{kprof}=PTi{kprof}*NaN;
    
    if ~isempty(Pt) || ~isempty(Ps),
        P = unique([Pt;Ps]);
        hi(kprof,7) = max(P);
        hi(kprof,8) = length(P);  
    end
    
end


% find bad location
I=find(hi(:,5)~=0);
hi=hi(I,:);
PTi=PTi(I); PSi=PSi(I);
PFi=PFi(I); POi=POi(I);
PLi=PLi(I);
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
PFi=PFi(I); POi=POi(I);
PLi=PLi(I);
hs=hs(I);

%save fcell
name_fcell=[conf.temporary_fcell info_deployment.EXP '_lr0_fcell.mat'];
save(name_fcell,'hi','hs','PTi','PSi','PFi','POi','PLi','EXP','PI','NATION','isoxy','isfluo','islight');

%% save in Argo netcdf format
if length(hi)>0
    suffix = 'lr0_prof.nc';
    convert_fcell2ARGO(conf,EXP,name_fcell,'lr0_prof.nc',[],one_smru_name);
end

disp(sprintf('\t%d tags',length(ltag)));
disp(sprintf('\t%d profiles',N));

