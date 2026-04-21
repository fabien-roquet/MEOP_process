%% sc_load_odv_fcell_argo

EXP=info_deployment.EXP;
PI=info_deployment.PI;
NATION=info_deployment.NATION;

fid = fopen([conf.rawdir info_deployment.nomfic]);
tline = fgetl(fid);
tline = fgetl(fid);
fclose(fid);

% read header
clear F O C BD;
isfluo=0; isoxy=0; iscond=0;islight=0;
Fluo=strfind(tline,'Fluorescence');
Oxy =strfind(tline,'Oxygen');
Cond=strfind(tline,'Conductivity');
Light=strfind(tline,'Light');

% read data
if length(Cond)>0 & length(Fluo)>0 & length(Oxy)>0    
   [smru_names,numProf,date,lon,lat,P,T,S,F,O,C] = ...
        textread([conf.rawdir info_deployment.nomfic],'%s%d%*s%s%f%f%*d%f%f%f%f%f%f',...
        'delimiter',';','headerlines',2);
    isfluo=double(length(find(F~=999))~=0); isoxy=double(length(find(O~=999))~=0);
    isS=find(S~=999);
    if length(isS)==0
        S=sw_salt_from_cond(C,T,P);
    end
    F(F==999)=NaN;
    O(O==999)=NaN;
    L=T.*NaN;
elseif length(Cond)>0
    [smru_names,numProf,date,lon,lat,P,T,S] = ...
        textread([conf.rawdir info_deployment.nomfic],'%s%d%*s%s%f%f%*d%f%f%f%*f',...
        'delimiter',';','headerlines',2);
    F=T.*NaN;
    O=T.*NaN;
    L=T.*NaN;   
elseif length(Fluo)==0 & length(Oxy)==0    
    [smru_names,numProf,date,lon,lat,P,T,S] = ...
        textread([conf.rawdir info_deployment.nomfic],'%s%d%*s%s%f%f%*d%f%f%f',...
        'delimiter',';','headerlines',2);
    F=T.*NaN;
    O=T.*NaN;
    L=T.*NaN;
elseif length(Fluo)>0 & length(Oxy)>0 & length(Light)>0
    [smru_names,numProf,date,lon,lat,P,T,S,F,L,O] = ...
        textread([conf.rawdir info_deployment.nomfic],'%s%d%*s%s%f%f%*d%f%f%f%f%f%f',...
        'delimiter',';','headerlines',2);
    isfluo=double(length(find(F~=999))~=0); isoxy=double(length(find(O~=999))~=0);islight=double(length(find(L~=999))~=0);
elseif length(Fluo)>0 & length(Oxy)>0
    [smru_names,numProf,date,lon,lat,P,T,S,F,O] = ...
        textread([conf.rawdir info_deployment.nomfic],'%s%d%*s%s%f%f%*d%f%f%f%f%f',...
        'delimiter',';','headerlines',2);
    isfluo=double(length(find(F~=999))~=0); isoxy=double(length(find(O~=999))~=0);
    L=T.*NaN;
elseif length(Light)>0
    [smru_names,numProf,date,lon,lat,P,T,S,F,L] = ...
        textread([conf.rawdir info_deployment.nomfic],'%s%d%*s%s%f%f%*d%f%f%f%f%f',...
        'delimiter',';','headerlines',2);
    isfluo=double(length(find(F~=999))~=0); islight=double(length(find(L~=999))~=0);
    O=T.*NaN;
end


% create header hi (fcell format)
Ihead=find(numProf~=0);
N=length(Ihead);
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

% split tags according to table_split_tags,
for ktag=1:length(ltag),
    smru_name = ltag{ktag};
    if any(strcmp(conf.table_split_tags.Properties.RowNames,smru_name)),
        nsplit = conf.table_split_tags{smru_name,'nsplit'};
        Iind = find(strcmp(hs,smru_name));
        jul = hi(Iind,4);
        [julsort,ii]=sort(diff(jul),'descend');
        Ibeg = [1;ii(1:nsplit-1)+1];
        Iend = [ii(1:nsplit-1);length(jul)];
        for ksplit = 1:nsplit,
            for kk = Ibeg(ksplit):Iend(ksplit)
                new_smru_name = sprintf('%s-N%1d',smru_name,ksplit);
                hs{Iind(kk)}=new_smru_name;
            end
        end
    end
end
ltag=unique(hs);
for ii=1:length(ltag),
    hi(strcmp(ltag{ii},hs),2)=ii;
end
hi(:,10)=hi(:,2);

% format data
PTi=cell(1,N); PSi=cell(1,N); 
PFi=cell(1,N); POi=cell(1,N);
PLi=cell(1,N);
T(T==999)=NaN; S(S==999)=NaN;
F(F==999)=NaN; O(O==999)=NaN;
L(L==999)=NaN;
Ihead2=[Ihead;length(P)+1];
nprof=1; ntag=1; pold=0;
for ii=1:N,    
    I=Ihead2(ii):Ihead2(ii+1)-1;
    I=I(diff([P(I);P(I(end))+1])~=0);
    PTi{ii}=[P(I) T(I)];
    PSi{ii}=[P(I) S(I)];
    PFi{ii}=[P(I) F(I)];
    POi{ii}=[P(I) O(I)];
    PLi{ii}=[P(I) L(I)];
    hi(ii,7)=P(Ihead2(ii+1)-1);
    hi(ii,8)=length(I);  
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
    convert_fcell2ARGO(conf,info_deployment.EXP,name_fcell,suffix,[],one_smru_name);
end

disp(sprintf('\t%d tags',length(ltag)));
disp(sprintf('\t%d profiles',N));

