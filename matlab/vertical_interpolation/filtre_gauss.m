function [PXi,Vi]=filtre_gauss(PX,ts,ind_p,Pi)
%     function [PXi,Vi]=filtre_gauss(PX,ts,ind_p,Pi)
%         
%     entree:  PX données à traiter; Pas de NaN dans PX!!
%       ts fenetre de lissage
%       ind_p indice de la colonne de référence dans PX (default =1)
%       Pi grille d'interpolation (defaut: PX(:,ind_p))
% 
%     sortie:  PXi signal interpolé aux points Pi.

%          Derniere modification: 18/02/07 Fabien Roquet

if nargin==3,
    Pi=PX(:,ind_p);
elseif nargin==2,
    ind_p=1;
    Pi=PX(:,ind_p);
end
Pi=reshape(Pi,length(Pi),1);

% classement par ordre croissant de pressions
P=PX(:,ind_p);
[P2,I]=sort(P);
PX2=PX(I,:);
[N1,N2]=size(PX2);
P2=PX2(:,ind_p);
PXi=zeros(length(Pi),N2);

for kk=1:length(Pi),
    
    dP=P2-Pi(kk);
    I=find(abs(dP)<=5*ts);
    if isempty(I) | all(dP(I)>0) | all(dP(I)<0)
        PXi(kk,:)=NaN;
    else
        f=exp((-dP(I).^2)/(1.44*ts.^2));
        d=PX2(I,:);
        PXi(kk,:)=sum((f*ones(1,N2)).*d,1)/sum(f);
    end
    
end

PXi(:,ind_p)=Pi;
Vi=Pi*0+1; %obsolete
