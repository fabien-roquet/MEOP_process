function[depthmax,vert_speed,da_duration,da_speed,delim]=poly_4_delim(tdr,chg,dureed) %#ok<*INUSD>
%% Délimitation du bottom des plongées
%
% On exploite la signature en sigmoïde du profil de la vitesse verticale.
% Pour cela on simplifie le signal avec un polynôme de degré 4. Les
% solutions de l'inégalité (I): -0.7 < P(x) < +0.7 correspondent aux lignes
% situées dans la phase de bottom. Il y a toutefois un bémol: si les
% limites du bottom trouvées de cette façon sont trop proches de la surface
% (ie. vérifie le critère suivant: profondeur < 40% profondeur maximum de
% la plongée) alors on effectue une correction. La nouvelle limite est
% alors la première ligne de prof > 40% de Pmax. Autre cas particulier, si
% aucune solution n'est trouvée à (I) alors on place les limites de part et
% d'autre de Pmax (ie. bottom presque inexistant).
%
display('    poly_4_delim... (3 steps)');
plotting=0;
resolution=round((tdr(2,1)-tdr(1,1))/(1.1574e-05));
displaytime=1;
mobmean=12;
crit=0.75;
%% 1°) Vitesse verticale
%
schg=size(chg);
s=size(tdr);
vert_speed=zeros(s(1),1);
for k=2:s(1)% vitesse vert du premier point = 0
    vert_speed(k)=(tdr(k,2)-tdr(k-1,2))/resolution;% vitesse en m/sec
end
vert_speed2=[];
vert_speed2(1:7)=vert_speed(1:7);
vert_speed2(end-6:end)=vert_speed(end-6:end);
for k=2+mobmean/2:s(1)-mobmean/2% lissage avec moyenne mobile pour
    %repérage des zones où la vitesse verticale est homogène
    vert_speed2(k)=mean(vert_speed(k-mobmean/2:k+mobmean/2));
end
vert_speed=vert_speed2';
clear vert_speed2;
display('    1°) done')
%% 2°) Extraction des profondeurs maximales
%
depthmax=zeros(2,schg(2));
for k=1:schg(2)
    [depthmax(1,k),ind]=max(tdr(chg(1,k):chg(2,k),2));
    depthmax(2,k)=ind+chg(1,k)-1;
end
display('    2°) done')
%% 3°) Identification de montées et descentes
%
delim=zeros(2,size(chg,2));
da_duration=zeros(2,size(chg,2));
da_speed=zeros(2,size(chg,2));
%
% initialisation du chronomètre
%
completed=1;
t=[];
tic
%
for s=1:size(chg,2)
    plg=vert_speed(chg(1,s):chg(2,s));
    %
    % fitting du polynome de degré 4
    %
    x_1 = (1:numel(plg))';% vecteur de la longueur de plg
    ok_ = isfinite(x_1) & isfinite(plg);
    ft_ = fittype('poly4');
    cf_ = fit(x_1(ok_),plg(ok_),ft_);% fitting de la vitesse verticale
    % en fonction du rang des valeurs.
    vf_ = feval(cf_,x_1);% calcul des valeurs prédites
    %
    % attribution des phases
    %
    ph=vf_>=crit;% -0.5 = descente = "valeurs prédites >= 0.4"
    ph=ph*-0.5;
    trans=vf_<-crit;% +0.5 = montée = "valeurs prédites <= -0.4"
    trans=trans*0.5;
    ph=ph+trans;% descentes et montées sont répertoriées dans "ph".
    % les éléments nuls de "ph" correspondent au bottom.
    clear trans
    ph(ph==0)=-1;% -1 = bottom = "-0.4 < valeurs prédites < 0.4"
    %
    % Cas particuliers
    %
    B=find(ph==-1);% une phase de bottom a t-elle été affectée ?
    A=find(ph==0.5);% une phase de montée a t-elle été affectée ?
    %% Correction résultat
    % 
    if isempty(B)==0% Il existe un bottom
        % verification que "prof debut du bottom > 80% prof max"
        %
        if tdr(B(1)+chg(1,s)-1,2)/depthmax(1,s)>=0.4% tout est en ordre
            delim(1,s)=B(1)+chg(1,s)-1;
        else% nouveau debut de bottom à 80% Pmax
            k=2;
            while tdr(B(k)+chg(1,s),2)/depthmax(1,s)<0.4 && k<length(B)
                k=k+1;
            end
            delim(1,s)=B(k)+chg(1,s)-1;
        end
    elseif isempty(B)==1% Il n'existe pas de bottom
        delim(1,s)=depthmax(2,s)-1;
    end
    if isempty(B)==0% la correction précédente a échouée
        if k>length(B)
            delim(1,s)=depthmax(2,s)-1;
        end
    end
    if isempty(A)==0% Il existe une montée
        % verification que "prof fin du bottom > 80% prof max"
        %
        if tdr(A(1)+chg(1,s)-1,2)/depthmax(1,s)>=0.4% tout est en ordre
            delim(2,s)=A(1)+chg(1,s)-1;
        else% nouvelle fin de bottom à 80% Pmax
            k=1;
            while tdr(A(1)-k+chg(1,s)-1,2)/depthmax(1,s)<0.4 && A(1)-k+chg(1,s)>delim(1,s)
                k=k+1;
            end
            delim(2,s)=A(1)-k+chg(1,s)-1;
        end
    elseif isempty(A)==1% Il n'existe pas de montée
        delim(2,s)=depthmax(2,s)+1;
    end
    if isempty(A)==0% la correction précédente a échouée
        if A(1)-k==delim(1,s)
            delim(2,s)=depthmax(2,s)+1;
        end
    end
    clear('B','A')
    %
    % Désormais "delim" stocke les statuts corrigés des phases de plongée
    %
    %% attribution du statut aux données tdr
    %
    tdr(chg(1,s):delim(1,s)-1,5)=-0.5;% -0.5 = descente
    tdr(delim(1,s):delim(2,s),5)=-1;% -1 = bottom
    tdr(delim(2,s)+1:chg(2,s),5)=0.5;% +0.5 = montée
    %
    tdr(chg(1,s):chg(2,s),6)=s;% ... et attribution du numéro de plongée
    da_duration(1,s)=(delim(1,s)-chg(1,s)+1)*resolution;%durée descente en sec
    da_duration(2,s)=(chg(2,s)-delim(2,s)+1)*resolution;%durée montée en sec
    da_speed(1,s)=mean(vert_speed(chg(1,s):delim(1,s)));%vitesse moy desc
    da_speed(2,s)=mean(vert_speed(delim(2,s):chg(2,s)));%vitesse moy asc
    if displaytime==1
        if s/size(chg,2)>0.1*completed
            t(completed)=toc; %#ok<AGROW>
            remaining_time=(10-completed)*t(completed)/completed;
            if completed==1
                current_time=clock;
                expected_end=disptime(current_time(4)*3600+current_time(5)*60+current_time(6)+remaining_time);
                remaining_time=disptime(remaining_time);
                str=sprintf('      remaining time: %d:%d:%d  ;  expected end at %d:%d:%d ',remaining_time(2),...
                    remaining_time(3),remaining_time(4),expected_end(2),expected_end(3),expected_end(4));
                disp(str)
            else
                remaining_time=disptime(remaining_time);
                str=sprintf('      remaining time: %d:%d:%d',remaining_time(2),...
                    remaining_time(3),remaining_time(4));
                disp(str)
            end
            completed=completed+1;
        end
    end
end
%
% graphique présentant le résultat
%
if plotting == 1
    figure('poly_4_delim result plot')
    plot(tdr(:,1),tdr(:,2)*-1,'k')
    hold all
    for s=1:size(chg,2)
        plot(tdr(delim(1,s),1),tdr(delim(1,s),2)*-1,'.g')
        plot(tdr(delim(2,s),1),tdr(delim(2,s),2)*-1,'.r')
    end
end
display('    3°) done')