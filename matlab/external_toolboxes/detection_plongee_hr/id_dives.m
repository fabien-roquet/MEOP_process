function[chg,tdr,y2,last,dureed]=id_dives(tdr,pseuil,plotting1,dureed)
%% Routine d'identification des plongées
%
% Attribution des statuts de surface et de bottom et supression des 
% périodes de bottom anormalement courtes.
%
%disp('    id_dives... (3 steps)')
%% 1°) Attribution des statuts "s"=1 et "b"=-1
%
resolution=round((tdr(2,1)-tdr(1,1))/(1.1574e-05)); %détermination de la 
% résolution des données en datenum, 1.1574e-05 correspond à 1 seconde
s=size(tdr);
states=zeros(s(1),1);
for k=1:s(1)
	if tdr(k,2)>=pseuil
		states(k)=-1;% --> bottom
	else
		states(k)=1;% --> surface
	end
end
divenumb=zeros(s(1),1);
tdr=cat(2,tdr,states,divenumb);% ajout des colonnes "states" et 
% "dive number"
clear('divenumb');
s=size(tdr); % mise a jour de "s"
%disp('    1°) done');
%% 2°) Identification des transition S/B et B/S
%
chg=[];
% chgtxt={'row numb','ptime','state'}; % Pour info
i=0;
t=0;% permet de commencer la séquence par une plongée et d'alterner avec 
% remontée obligatoirement
for k=1:s(1)-1
    if states(k)-states(k+1)==2 && t==0
        i=i+1;
        t=1;
        chg(1,i)=k; %#ok<AGROW> % plongées en premiere ligne
    elseif states(k)-states(k+1)==-2 && t==1
        t=0;
        chg(2,i)=k; %#ok<AGROW> % remontées en deuxieme ligne
    end	
end
clear('states');
schg=size(chg);
last=chg(2,schg(2)-1);
t=0;
if chg(2,schg(2))==0% supression de la dernière ligne si pas d'info sur la 
    %remontée
    last=chg(1,schg(2));% mais enristrement de cette valeur car on peut
    % savoir le temps passé en surface après la dernière plongée
    t=1;% dans ce cas le tes du 3°) porte aussi sur cette valeur
    chg(:,schg(2))=[];
    schg=size(chg);% mise a jour de "schg"
end
%disp('    2°) done');
%% 3°) Suppression des périodes de plg trop courtes (< duree)
%
if plotting1==1
    shortdives=[];
    i=1;
    for k=1:schg(2)
        if (chg(2,k)-chg(1,k))*resolution<=500
            shortdives(i)=(chg(2,k)-chg(1,k))*resolution; %#ok<AGROW>
            i=i+1;
        end
    end
    figure1 = figure('PaperType','a4letter','PaperSize',[20.98 29.68],...
    'Name','ana_dives: plot1.2');
    axes('Parent',figure1,'YGrid','on',...
    'XTickLabel',...
{'0-50','50-100','100-150','150-200','200-250','250-300','300-350','350-400','400-450','450-500'},...
    'XTick',[0 50 100 150 200 250 300 350 400 450 500]);
    xlim([-27.8 527.8]);
    ylim([0 10]);
    box('on');
    hold('all');
    xlabel('dive time in seconds');
    ylabel('number');
    abscisse=0:50:500;
    histo=histc(shortdives,abscisse);
    bar(abscisse,histo)
    title 'Short dive periods histogram'
%disp('Rentrez une valeur minimale de temps (en secondes) en dessous de la quelle');
dureed=input('les périodes de plongée seront supprimées (conseillé 500): ');
    clear('shortdives','histo','abscisse');
end
y2=0;
k=1;
while k<=schg(2);
    if k>schg(2)
        break
    end
    if ((chg(2,k)-chg(1,k))*resolution)<dureed
        for i=chg(1,k):chg(2,k)% correction dans "tdr"
            tdr(i,s(2)-1)=1;
        end
        % suppression de la plongée dans chg
        chg(:,k)=[]; %#ok<AGROW>
        schg=size(chg);
        y2=y2+1;% juste pour infos dans 7°)
    else
        k=k+1;
    end
end
%disp('    3°) done');