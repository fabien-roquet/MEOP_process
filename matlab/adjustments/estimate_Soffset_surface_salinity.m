function Soffset = estimate_Soffset_surface_salinity(argo_qc,zone)
% compare surface salinity 
% By default, comparison is made along the Antarctic continental shelf (zone='ACS')
%             zone='NKP': comparison over the Northern Kerguelen Plateau

if nargin==1,
    zone='PB'; % Prydz Bay
end

switch zone
    case 'PB' % Prydz Bay
        imask = 5;
        lambdax=0.5;lambday=.25;lambdat=15;
        lat_bin=-70:.02:-66; lon_bin=65:.04:85;
        lat_bin_centre=lat_bin+.01;
        lon_bin_centre=lon_bin+.02;
    case 'NKP' % Northern Kerguelen Plateau
        imask = 6;
        lambdax=1;lambday=.5;lambdat=15;
        lat_bin=-55:.04:-50; lon_bin=70:.08:78;
        lat_bin_centre=lat_bin+.02;
        lon_bin_centre=lon_bin+.04;
    case 'CS' % Casey Station
        imask = 6;
        lambdax=.5;lambday=.25;lambdat=15;
        lat_bin=-67:.02:-65.2; lon_bin=105:.04:111;
        lat_bin_centre=lat_bin+.01;
        lon_bin_centre=lon_bin+.02;
end

load('mask.mat')
argo_qc.mask = interp2(lonmask,latmask,mask,argo_qc.LONGITUDE,argo_qc.LATITUDE,'nearest');
if strcmp(zone,'CS'), argo_qc.mask=6; end

Nm=[]; Xm=[]; Ym=[]; Tm=[]; Sm=[];
for ii=1:argo_qc.ntag,
    
    descr=argo_qc.list_descr{ii};
    K=find(strcmp(argo_qc.platform_number,descr) & argo_qc.mask==imask & ~isnan(argo_qc.PSAL(2,:)'));
    if length(K)>20,
        Nm=[Nm; ii*ones(length(K),1)    ];
        Xm=[Xm; argo_qc.LONGITUDE(K)    ];
        Ym=[Ym; argo_qc.LATITUDE(K)     ];
        Tm=[Tm; argo_qc.JULD_LOCATION(K)];
        Sm=[Sm; argo_qc.PSAL(2,K)'      ];
    end

end

time_bin=min(Tm):4:max(Tm)+4;
time_bin_centre=time_bin+2;

liste_tag=unique(Nm);
disp(sprintf('%d profils mirounga utilisés',length(Nm)));


%% Rééchantillonnage des Sz sur une grille régulière

count=0;
Sm2=[];Nm2=[];Xm2=[];Ym2=[];Tm2=[];

for nseal=1:length(liste_tag),
    
    K=find(Nm==liste_tag(nseal));
    [n,bin] = histc(Ym(K),lat_bin);
    for kk=1:length(n)-1,
        if n(kk)==0, continue, end
        I=find(bin==kk);
        [n2,bin2] = histc(Xm(K(I)),lon_bin);
        for kk2=1:length(n2)-1
            if n2(kk2)==0, continue, end
            J=find(bin2==kk2);
            [n3,bin3] = histc(Tm(K(I(J))),time_bin);
            for kk3=1:length(n3)-1
                if n3(kk3)==0, continue, end
                L=find(bin3==kk3);
                count=count+1;
                Nm2(count)=liste_tag(nseal);
                Xm2(count)=lon_bin_centre(kk2);
                Ym2(count)=lat_bin_centre(kk);
                Tm2(count)=time_bin_centre(kk3);
                Sm2(count)=mean(Sm(K(I(J(L)))));
            end
        end
    end
end

Sm=Sm2';
Nm=Nm2';
Xm=Xm2';
Ym=Ym2';
Tm=Tm2';
liste_tag=unique(Nm);
disp(sprintf('%d obs mirounga after re-sampling',length(Sm2)));

%% Intercomparaison entre les donnees S0 mir07
data=Sm; xdata=Xm; ydata=Ym; tdata=Tm; ndata=Nm;

% signal variance
varsig=0;
for nseal=1:length(liste_tag),
    varsig=varsig+var(data(ndata==liste_tag(nseal)));
end
varsig=varsig/length(liste_tag);
varnoise=0.005^2;

% covariance des obs
Cdata=-(repmat(xdata,1,length(xdata))-repmat(xdata,1,length(xdata))').^2/lambdax^2;
Cdata=Cdata-(repmat(ydata,1,length(ydata))-repmat(ydata,1,length(ydata))').^2/lambday^2;
Cdata=Cdata-(repmat(tdata,1,length(tdata))-repmat(tdata,1,length(tdata))').^2/lambdat^2;
Cdata=varsig*exp(Cdata)+varnoise*eye(length(Cdata));

% retrait des observations isolées
Iiso=[];
for kk=1:length(ndata),
    I=find(ndata~=ndata(kk));
    Iiso_=setdiff(1:length(ndata),Iiso);
    if all(Cdata(I,kk)<.3*varsig) & length(find(ndata(Iiso_)==ndata(kk)))>10, 
        Iiso(end+1)=kk; 
    end
end
data(Iiso)=[];xdata(Iiso)=[];ydata(Iiso)=[];tdata(Iiso)=[];ndata(Iiso)=[];
liste_tag_i=[];
for nseal=1:length(liste_tag),
    if ~isempty(find(ndata==liste_tag(nseal))),
        liste_tag_i(end+1)=liste_tag(nseal);
    end
end
Cdata(Iiso,:)=[]; Cdata(:,Iiso)=[];
disp(sprintf('%d obs mirounga after removal of isolated observations',length(data)));

%% Winv weighting matrix
n=length(data)-1;  
Winv=zeros(n);
I=zeros(n,1); J=I; BI=I; BJ=I;

for ii=1:n,
    I(ii)=ii+1;
    J(ii)=1;
    BI(ii)=ndata(I(ii));
    BJ(ii)=ndata(J(ii));
    Winv(ii,1:ii)= Cdata(I(ii),I(1:ii)) + Cdata(J(ii),J(1:ii)) ...
        - Cdata(I(ii),J(1:ii)) - Cdata(J(ii),I(1:ii)) ;
    Winv(1:ii,ii)= Winv(ii,1:ii)';
end

Y=zeros(n,1); dlon=Y; dlat=Y; dtime=Y;
M=length(liste_tag_i); X=zeros(n,M);
for nn=1:n,
    Y(nn)=data(J(nn))-data(I(nn));
    dlon(nn)=xdata(J(nn))-xdata(I(nn));
    dlat(nn)=ydata(J(nn))-ydata(I(nn));
    dtime(nn)=tdata(J(nn))-tdata(I(nn));
end
for kk=1:length(liste_tag_i),
    for nn=1:n,
        if BI(nn)==liste_tag_i(kk), X(nn,kk)=1; end
        if BJ(nn)==liste_tag_i(kk), X(nn,kk)=-1; end
    end
end

% choix des balises prises en compte dans l'analyse
liste_tag_f=liste_tag_i;
X2=X(:,2:end);

% Calcul de la solution
[L U]=lu(Winv);
A=inv((X2'/U)*(L\X2));
beta=A*(X2'/U)*(L\Y);

% Calcul du vecteur de données corrigées
Sol=[0;beta];cdata=data*0;
for kk=1:length(liste_tag_f),
    cdata(ndata==liste_tag_f(kk))=data(ndata==liste_tag_f(kk))+Sol(kk);
end

% % affichage des statistiques
% disp(' ');
% for kk=1:length(liste_tag_f),
%     I=find(ndata==liste_tag_f(kk));
%     ds2=[];
%     for ii1=1:length(I),
%         C=exp( (xdata(I)-xdata(I(ii1))).^2/lambdax^2 + ...
%                (ydata(I)-ydata(I(ii1))).^2/lambday^2 + ...
%                (tdata(I)-tdata(I(ii1))).^2/lambdat^2  );
%         C(ii1)=NaN;
%         [c,ii2]=min(C);
%         ds2(ii1)=(data(I(ii2))-data(I(ii1)))^2;
%     end
%     disp( sprintf( '%2d %3d %4.3f %4.3f %4.3f %4.3f %4.3f', ...
%         liste_tag_f(kk), ...
%         length(I), ...
%         mean( data(I) ), ...
%         std( data(I) ), ...
%         sqrt(nanmean(ds2)/2) , ...
%         Sol(kk), ...
%         mean( cdata(I) )));
% end
% disp(sprintf('mean %5.3f ; std %5.3f',mean(Sol),std(Sol)))
% 
% vardiff=[];
% vardiff(2:length(A)+1,2:length(A)+1) = ...
%     ones(length(A),1)*diag(A)'+diag(A)*ones(1,length(A),1)-2*A;
% vardiff(2:length(A)+1,1)=diag(A);
% vardiff(1,2:length(A)+1)=diag(A)';
% stddiff=sqrt(vardiff);
% disp(' ');
% fmt=''; for kk=1:length(liste_tag_f), fmt(end+1:end+6)='%4.0f '; end;
% disp( sprintf(fmt,liste_tag_f) );
% fmt=''; for kk=1:length(liste_tag_f), fmt(end+1:end+6)='%4.1f '; end; fmt=[fmt '\n'];
% disp( sprintf(fmt,stddiff*1000) );

%% output result

Soffset = zeros(argo_qc.ntag,1);
Soffset(liste_tag_f) = Sol - mean(Sol);

for ii=1:argo_qc.ntag,
        
    descr=argo_qc.list_descr{ii};
    disp(sprintf('%5.2f', Soffset(ii)));
    
end


