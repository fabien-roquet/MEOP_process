%% ct3
load('/Users/roquet/GoogleDrive/MEOP-CTD/original_datasets/ct3_cor.mat')
I=find(ismember(hi(:,10),19+[1:10]));

hi=hi(I,:);
PTi=PTi(I);
PSi=PScor(I);
PFi=PTi; for kk=1:length(PFi), PFi{kk}=[]; end
POi=PFi;
EXP='ct3';
PI='XTOPHE';
isoxy=0;
isfluo=0;

[ltag,ia,ic] = unique(hi(:,10));
hi(:,10)= ic;
hs={};
seal_name= {
    'ct3-9927-04'
    'ct3-9935-04'
    'ct3-9938-04'
    'ct3-9931-04'
    'ct3-9932-04'
    'ct3-9934-04'
    'ct3-9928-04'
    'ct3-9929-04'
    'ct3-9926-04'
    'ct3-9118-04'
    };
for kk=1:length(ic),
    hs{kk}=seal_name{ic(kk)};
end
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
hs=hs(I);

name_fcell=[conf.rawdir EXP '_fcell.mat'];
save(name_fcell,'hi','hs','PTi','PSi','PFi','POi','EXP','PI','isoxy','isfluo');

%% ct7
load('/Users/roquet/GoogleDrive/MEOP-CTD/original_datasets/ct7_ct11_cor.mat')
I=find(ismember(hi(:,10),29+[1:5 10:11]));

hi=hi(I,:);
PTi=PTi(I);
PSi=PScor(I);
for kk=1:length(PSi),
    if ~isempty(PSi{kk}),
        PSi{kk}(:,3)=[];
    end
end
PFi=PTi; for kk=1:length(PFi), PFi{kk}=[]; end
POi=PFi;
EXP='ct7';
PI='XTOPHE';
isoxy=0;
isfluo=0;

[ltag,ia,ic] = unique(hi(:,10));
hi(:,10)= ic;
hs={};
seal_name= {
    'ct7-10028-05'
    'ct7-10029-05'
    'ct7-10037-05'
    'ct7-10035-05'
    'ct7-10026-05'
    'ct7-10032-05'
    'ct7-10039-05'
    };
for kk=1:length(ic),
    hs{kk}=seal_name{ic(kk)};
end
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
hs=hs(I);

name_fcell=[conf.rawdir EXP '_fcell.mat'];
save(name_fcell,'hi','hs','PTi','PSi','PFi','POi','EXP','PI','isoxy','isfluo');

%% ct11
load('/Users/roquet/GoogleDrive/MEOP-CTD/original_datasets/ct7_ct11_cor.mat')
I=find(ismember(hi(:,10),[35 37 38]));

hi=hi(I,:);
PTi=PTi(I);
PSi=PScor(I);
for kk=1:length(PSi),
    if ~isempty(PSi{kk}),
        PSi{kk}(:,3)=[];
    end
end
PFi=PTi; for kk=1:length(PFi), PFi{kk}=[]; end
POi=PFi;
EXP='ct11';
PI='XTOPHE';
isoxy=0;
isfluo=0;

[ltag,ia,ic] = unique(hi(:,10));
hi(:,10)= ic;
hs={};
seal_name= {
    'ct11-10064-05'
    'ct11-10096-05'
    'ct11-10097-05'
    };
for kk=1:length(ic),
    hs{kk}=seal_name{ic(kk)};
end
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
hs=hs(I);

name_fcell=[conf.rawdir EXP '_fcell.mat'];
save(name_fcell,'hi','hs','PTi','PSi','PFi','POi','EXP','PI','isoxy','isfluo');

%% wd3
load('/Users/roquet/GoogleDrive/MEOP-CTD/original_datasets/wd3_cor.mat')
I=find(ismember(hi(:,10),1:3));

hi=hi(I,:);
PTi=PTi(I);
PSi=PSi(I);
PFi=PTi; for kk=1:length(PFi), PFi{kk}=[]; end
POi=PFi;
EXP='wd3';
PI='HINDELL';
isoxy=0;
isfluo=0;

[ltag,ia,ic] = unique(hi(:,10));
hi(:,10)= ic;
hs={};
seal_name= {
    'wd3-CTD1-07'
    'wd3-CTD2-07'
    'wd3-CTD3-07'
    };
for kk=1:length(ic),
    hs{kk}=seal_name{ic(kk)};
end
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
hs=hs(I);

name_fcell=[conf.rawdir EXP '_fcell.mat'];
save(name_fcell,'hi','hs','PTi','PSi','PFi','POi','EXP','PI','isoxy','isfluo');

