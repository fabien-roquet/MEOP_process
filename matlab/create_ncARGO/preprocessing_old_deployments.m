%% for early deployments with two ODV version, concatenate files keeping only highest resolution data
list_ex={'ct10','ct12','ct17'};
for ii=1:length(list_ex)
    if ~exist([conf.rawdir list_ex{1,ii} 'r' '_ODV.txt'],'file'),
        continue,
    end
    fid = fopen([conf.rawdir list_ex{1,ii} 'r' '_ODV.txt']);
    tline = fgetl(fid);
    tline = fgetl(fid);
    fclose(fid);
    J=find(strcmp([list_ex{1,ii}],list_deployment{:,1}));
    [tagr,numProfr,typer,dater,lonr,latr,Dr,Pr,Tr,Sr,Fr,Or,Cr] = ...
        textread([conf.rawdir list_ex{1,ii} 'r' '_ODV.txt'],'%s%s%s%s%s%s%s%f%f%f%f%f%f',...
        'delimiter',';','headerlines',2);
    copyfile([conf.rawdir list_ex{1,ii} 'r' '_ODV.txt'],[conf.rawdir list_ex{1,ii} 'r' 'save_ODV.txt'],'f');
    [tag,numProf,type,date,lon,lat,D,P,T,S,F,O,C] = ...
        textread([conf.rawdir list_deployment{1,1}{J} '_ODV.txt'],'%s%s%s%s%s%s%s%f%f%f%f%f%f',...
        'delimiter',';','headerlines',2);
    copyfile([conf.rawdir list_deployment{1,1}{J} '_ODV.txt'],[conf.rawdir list_deployment{1,1}{J} 'save_ODV.txt'],'f');
    tagr=regexprep(tagr,'r','');
    indr=unique(tagr(~strcmp(tagr,'')));
    ind=unique(tag(~strcmp(tag,'')));
    U=intersect(indr,ind);
    fileID = fopen([conf.rawdir list_ex{1,ii} 'temp_ODV.txt'],'w');
    tline2='Cruise;Station;Type;yyyy-mm-dd hh:mm;Longitude [degrees_east];Latitude [degrees_north];Bot. Depth [m];Pressure [dbar];Temperature [C];Salinity [PSU];Fluorescence;Oxygen [umol/l]';
    fprintf(fileID,'%s\n','// Generic ODV file');
    fprintf(fileID,'%s\n',tline2);
    Tr2=Tr; Tr2(Tr2<-3|Tr2>40)=NaN;
    Cr2=Cr; Cr2(Cr2<10|Cr2>50)=NaN;
    Sr=sw_salt_from_cond(Cr2,Tr2,Pr); Sr(isnan(Sr))=999;
    for jj=1:length(U)
        I=find(strcmp(U(jj),tagr));
        fin=I(end)+1;
        while fin <=length(Tr) & strcmp(tagr(fin),'')
            fin=fin+1;
        end
        for kk=I(1):fin-1
            fprintf(fileID,'%s;%s;%s;%s;%s;%s;%s;%3.4f;%3.4f;%3.4f;%3.4f;%3.4f\n',...
                tagr{kk},numProfr{kk},typer{kk},dater{kk}...
                ,lonr{kk},latr{kk},Dr{kk},Pr(kk),Tr(kk)...
                ,Sr(kk),Fr(kk),Or(kk));
        end
    end
    U=setdiff(ind,indr);
    T2=T; T2(T2<-3|T2>40)=NaN;
    C2=C; C2(C2<10|C2>50)=NaN;
    S=sw_salt_from_cond(C2,T2,P); S(isnan(S))=999;
    for jj=1:length(U)
        I=find(strcmp(U(jj),tag));
        fin=I(end)+1;
        while fin <=length(T) & strcmp(tag(fin),'')
            fin=fin+1;
        end
        for kk=I(1):fin-1
            fprintf(fileID,'%s;%s;%s;%s;%s;%s;%s;%3.4f;%3.4f;%3.4f;%3.4f;%3.4f\n',...
                tag{kk},numProf{kk},type{kk},date{kk}...
                ,lon{kk},lat{kk},D{kk},P(kk),T(kk)...
                ,S(kk),F(kk),O(kk));
        end
    end
    fclose(fileID);
    copyfile([conf.rawdir list_deployment{1,1}{J} 'temp_ODV.txt'],[conf.rawdir list_deployment{1,1}{J} '_ODV.txt'],'f');
    delete([conf.rawdir list_deployment{1,1}{J} 'temp_ODV.txt']);
    delete([conf.rawdir list_deployment{1,1}{J} 'r_ODV.txt']);
    
end

%% for early deployment ct1, merge temperature and salinity profiles
list_ex={'ct1','ct6'};
for ii=1:length(list_ex)
    if ~exist([conf.rawdir list_ex{1,ii} '_ODV.txt'],'file') | exist([conf.rawdir list_ex{1,ii} 'save_ODV.txt'],'file'),
        continue,
    end
    fid = fopen([conf.rawdir list_ex{1,ii} '_ODV.txt']);
    tline = fgetl(fid);
    tline = fgetl(fid);
    fclose(fid);
    [tag,numProf,date,lon,lat,P,T,S] = ...
        textread([conf.rawdir list_ex{1,ii} '_ODV.txt'],'%s%s%*s%s%f%f%*s%f%f%f%*f%*f',...
        'delimiter',';','headerlines',2);
    copyfile([conf.rawdir list_ex{1,ii} '_ODV.txt'],[conf.rawdir list_ex{1,ii} 'save_ODV.txt'],'f');
    fileID = fopen([conf.rawdir list_ex{1,ii} '_ODV.txt'],'w');
    fprintf(fileID,'%s\n','// Generic ODV file');
    fprintf(fileID,'%s\n',tline);
    nheadline=[find(~strcmp(tag,''));length(tag)+1];
    jj=1; nprof=1;
    while jj<length(nheadline)-1
        if strcmp(date(nheadline(jj)),date(nheadline(jj+1))),
            I1=nheadline(jj):nheadline(jj+1)-1;
            I2=nheadline(jj+1):nheadline(jj+2)-1;
            [Psort1,Isort1]=sort(P(I1));
            [Psort2,Isort2]=sort(P(I2));
            [Pinters,K1,K2]=intersect(Psort1,Psort2);
            T(I1(K1))=T(I2(K2));
            Itot=setdiff([I1 I2],I2(K2));
            [Psort,Itot2]=sort(P(Itot));
            Itot=Itot(Itot2);
            fprintf(fileID,'%s;%d;B;%s;%3.4f;%3.4f;0;%3.4f;%3.4f;%3.4f;999.0000;999.0000\n',...
                tag{nheadline(jj)},nprof,date{nheadline(jj)},lon(nheadline(jj)),lat(nheadline(jj)),P(Itot(1)),T(Itot(1)),S(Itot(1)));
            for kk=2:length(Itot)
                fprintf(fileID,';;;;;;;%3.4f;%3.4f;%3.4f;999.0000;999.0000\n',P(Itot(kk)),T(Itot(kk)),S(Itot(kk)));
            end
            jj=jj+2; nprof=nprof+1;
        else
            Itot=nheadline(jj):nheadline(jj+1)-1;
            fprintf(fileID,'%s;%d;B;%s;%3.4f;%3.4f;0;%3.4f;%3.4f;%3.4f;999.0000;999.0000\n',...
                tag{nheadline(jj)},nprof,date{nheadline(jj)},lon(nheadline(jj)),lat(nheadline(jj)),P(Itot(1)),T(Itot(1)),S(Itot(1)));
            for kk=2:length(Itot)
                fprintf(fileID,';;;;;;;%3.4f;%3.4f;%3.4f;999.0000;999.0000\n',P(Itot(kk)),T(Itot(kk)),S(Itot(kk)));
            end
            jj=jj+1; nprof=nprof+1;
        end
    end
    fclose(fileID);
end

