function[sinu]=sinuosity(tdr,depthmax,vert_speed,daindexes)
%% calcul de la sinuosité verticale des plongées
%
% sinuosité = bottomvertdistance observed / bottomvertdistance euclidian
%
%disp('    sinuosity...')
sda=size(daindexes);
euclidian=zeros(1,sda(2));
for k=1:sda(2)
    euclidian(k)=(depthmax(1,k)-tdr(daindexes(1,k),2))+(depthmax(1,k)-tdr(daindexes(2,k),2));
end
observed=zeros(1,sda(2));
sinu=zeros(1,sda(2));
for k=1:sda(2)
    local_extr=[]; % vecteur contenant les rg/indices des extremum locaux 
    % de la profondeur du bottom en fonction du temps
    j=1;
    for i=daindexes(1,k):daindexes(2,k)-1
        % si la vitesse verticale change de signe entre le rg i le rg i+1
        % alors on considère un extremum local en i
        if vert_speed(i)*vert_speed(i+1)<0
            local_extr(j)=tdr(i,2); %#ok<AGROW>
            j=j+1;
        end
    end
    if isempty(local_extr)==0 % si il existe au moins 1 extremum
        % initialisation de obs: distance verticale entre premier extr. et
        % le debut du bottom
        % +
        % distance verticale entre dernier extr. et la
        % fin du bottom
        observed(k)=abs(local_extr(1)-tdr(daindexes(1,k),2))+abs(local_extr(length(local_extr))-tdr(daindexes(2,k),2));
    else
        % bottom parfaitement plat
        % impossible mais bon, il faut tjr prévoir le pire.
        sinu(k)=1;
    end
    if length(local_extr)>1
        for i=2:length(local_extr)
            % somme des distances (tjr positives !) verticales parurues
            observed(k)=observed(k)+abs(local_extr(i)-local_extr(i-1));
        end
    end
    sinu(k)=observed(k)/euclidian(k);
end

