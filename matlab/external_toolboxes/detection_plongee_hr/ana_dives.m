function [dives,info_ana_dives,varargout]=ana_dives(tdr)
%%  Extraction des plongees et calcul de leurs parametres
%
% [x,y,z,t,a,b,c,d,e,f]=dives(tdr)
%
% x est le tableau de synthese des plongees a l'issue de l'analyse des
% donnees TDR
% y est une structure contenant des informations sur l'identification des
% plongees.
% z, facultatif - noms des colonnes de x.
% t, facultatif - les donnees tdr avec le statut (surf/plg) mis a jour.
% a, facultatif - noms des colonnes de t.
% b, facultatif - indexes de debut et fin de plongee (tableau "chg").
% c, facultatif - indexes de debut et fin de bottom (tableau "daindexes").
% e, facultatif - vitesse verticale lissee.
% f, facultatif - profondeur max et indexes.
%
% *Parametres:*
%
% * pseuil: profondeur a partir de laquelle l'animal est considere en
% plongee.
%
% * startnumb: premier num pour numerotation des plongees.
%
% * plotting1: "1" pour tracer l'histogramme des petites periode de surface 
% et de plongee. % Il permet a l'utilisateur de voir quel est le parametre
% "duree" le plus adapte aux donnees dans chacun des cas.
%
% * dureed: dans le cas ou l'histogramme n'est pas affiche, la valeur par
% defaut du parametre "duree" pour les periodes de plongee.
% Il s'agit du minimum de temps (en seconde) que doit durer une periode de 
% plongee pour etre considere comme valide.
%
% * infos: "1" pour afficher les infos a la fin de la fonction.
%
pseuil=15;
startnumb=1;
plotting1=0;
dureed=300;
display('ana_dives... (4 steps)');
%% 1) Identification des phases des plongees
%
resolution=(tdr(2,1)-tdr(1,1))/(1.1574e-05); %determination de la 
% resolution des donnees en datenum, 1.1574e-05 correspond a 1 seconde
[chg,tdr,y2,last,dureed]=id_dives(tdr,pseuil,plotting1,dureed);
s=size(tdr); % mise a jour de "s"
tdrtxt={'Ptime Num','Corrected Depth','External Temperature',...
    'Light Level','State','Dive Number'};
varargout(3)={tdrtxt};
display('1) done');
%% 2) Identification de montees et descentes
%
[depthmax,vert_speed,daduration,daspeed,daindexes]=poly_4_delim(tdr,chg,dureed);
[sinu]=sinuosity(tdr,depthmax,vert_speed,daindexes);
display('2) done');
%% 3) Affectation des valeurs dans la matrice 'dives'
%
schg=size(chg);
dives=zeros(schg(2),13);
divestxt={'dive numb','start ptime','start index','dive duration','max depth',...
    'time spent in surface after diving',...
    'desc duration','bottom duration','asc duration','bottom duration residual',...
    'desc speed','asc speed','sinuosity'};
varargout(1)={divestxt};
for k=1:schg(2)
	dives(k,1)=startnumb-1+k; % dive numb
    for i=chg(1,k)+1:chg(2,k)
        tdr(i,s(2))=dives(k,1);
    end
	dives(k,2)=tdr(chg(1,k)+1,1); % start ptime
    dives(k,3)=chg(1,k)+1; % 'start index
	dives(k,4)=(chg(2,k)-chg(1,k))*resolution; % dive duration
	dives(k,5)=depthmax(1,k); % max depth
    if k<schg(2)
		dives(k,6)=(chg(1,k+1)-chg(2,k))*resolution; % time spent in surface after diving
    else
        dives(k,6)=(last-chg(2,k))*resolution; % idem mais pour la derniere plongee
    end
	dives(k,7)=daduration(1,k); % desc duration
	dives(k,8)=dives(k,4)-daduration(1,k)-daduration(2,k); % bottom duration
	dives(k,9)=daduration(2,k); % asc duration
    dives(k,11)=daspeed(1,k); % desc speed
    dives(k,12)=daspeed(2,k); % asc speed
    dives(k,13)=sinu(k); % sinuosity
end
%
% regression lineaire multiple pour calcul des residus du bottom time
% model: Bottom time ~ maxdepth + dive duration
%
[betahat,Ibeta,res]=regress(dives(:,8),dives(:,4:5)); %#ok<ASGLU>
for k=1:schg(2)
    dives(k,10)=res(k); % bottom duration residual
end
display('3) done');
%% 4) Information sur le deroulement de la fonction
%
% & ajouts d'arguments de sortie
%
varargout(2)={tdr};
varargout(4)={chg};
varargout(5)={daindexes};
varargout(6)={vert_speed};
varargout(7)={depthmax};
varargout(8)={resolution};
clear('vert_speed','chg','daduration','depthmax');
info_ana_dives=struct('Seuil_prof', pseuil,...
    'Resolution_TDR',floor(resolution),...
    'Duree_minimale_plongee',dureed,...
    'Nb_plongee_reecrites_en_surface',y2) %#ok<NOPRT>
display('4) done');

