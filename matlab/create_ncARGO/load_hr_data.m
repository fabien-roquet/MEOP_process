function hrdata = load_hr_data(name_hr_file,continuous)


    hrdata = [];
    
    %% chargement du fichier ctd haute resolution
    fid = fopen(name_hr_file);
    tline = fgetl(fid);
    fclose(fid);
    
    Fluo = [strfind(tline,'Fluo') strfind(tline,'FLUORO')];
    Oxy = strfind(tline,'Oxy');
    Light = strfind(tline,'PPFD');
    % sur les nouveaux fichiers hr il y a 3 colonnes de pressions en plus
    newfile = length(strfind(tline,'PRESSURE_RE_SURFACE'));
    hrdata.isfluo=length(Fluo); hrdata.isoxy=length(Oxy); hrdata.islight=length(Light);
    hrdata.continuous = continuous;
    
    try
        if ~hrdata.isfluo & ~hrdata.isoxy & ~hrdata.islight & ~newfile % ancien fichier HR ctd
            [str_date,P,T,S] = ...
                textread(name_hr_file,'%s%f%f%*f%f%*[^\n]',...
                'delimiter','\t','headerlines',1);
            hrdata.F=T.*NaN;
            hrdata.O=T.*NaN;
            hrdata.L=T.*NaN;
            
        elseif  ~hrdata.isfluo & ~hrdata.isoxy & ~hrdata.islight & newfile % nouveau fichier HR ctd
            [str_date,P,T,S] = ...
                textread(name_hr_file,'%s%f%*f%*f%f%*f%f%*[^\n]',...
                'delimiter','\t','headerlines',1);
            hrdata.F=T.*NaN;
            hrdata.O=T.*NaN;
            hrdata.L=T.*NaN;
        
        elseif hrdata.isfluo & ~hrdata.isoxy & ~hrdata.islight
            [str_date,P,T,S,F] = ...
                textread(name_hr_file,'%s%f%f%*f%f%f%*[^\n]',...
                'delimiter','\t','headerlines',1);
            hrdata.isfluo=double(length(find(F~=999))~=0);
            hrdata.F=F;
            hrdata.O=T.*NaN;
            hrdata.L=T.*NaN;
        
        elseif ~hrdata.isfluo & hrdata.isoxy & ~hrdata.islight
            [str_date,P,T,S,O] = ...
                textread(name_hr_file,'%s%f%f%*f%f%f%*[^\n]',...
                'delimiter','\t','headerlines',1);
            hrdata.isoxy=double(length(find(O~=999))~=0);
            hrdata.F=T.*NaN;
            hrdata.O=O;
            hrdata.L=T.*NaN;
            
        elseif strfind(name_hr_file,'2016/14331_ctd.txt') % ct132-331-16,14331,2016
            [str_date,P,T,S,L] = ...
                textread(name_hr_file,'%s%f%f%*f%f%f%*[^\n]',...
                'delimiter','\t','headerlines',1);
            hrdata.islight=1;
            hrdata.F=T.*NaN;
            hrdata.O=T.*NaN;
            hrdata.L=L;

        elseif strfind(name_hr_file,'2014/13034_ctd.txt') %ct112-034-14,13034,2014
            [str_date,P,T,S] = ...
                textread(name_hr_file,'%s%f%*f%*f%f%*f%f%*[^\n]',...
                'delimiter','\t','headerlines',1);
            hrdata.F=T.*NaN;
            hrdata.O=T.*NaN;
            hrdata.L=T.*NaN;
        
        elseif ~hrdata.isfluo & ~hrdata.isoxy & hrdata.islight 
            [str_date,P,T,S,L] = ...
                textread(name_hr_file,'%s%f%*f%*f%f%*f%f%f%*f%*f%*[^\n]',...
                'delimiter','\t','headerlines',1);
            hrdata.islight=1;
            hrdata.F=T.*NaN;
            hrdata.O=T.*NaN;
            hrdata.L=L;
            
        elseif hrdata.isfluo & ~hrdata.isoxy & hrdata.islight 
            [str_date,P,T,S,F,L] = ...
                textread(name_hr_file,'%s%f%*f%*f%f%*f%f%f%f%*f%*[^\n]',...
                'delimiter','\t','headerlines',1);
            hrdata.isfluo=1;
            hrdata.islight=1;
            hrdata.F=F;
            hrdata.O=T.*NaN;
            hrdata.L=L;
                  
        end
    catch
        disp(['Error reading file: ' name_hr_file]);
    end
    
    % correction des bugs temporelles des nouvelles CTD, dans le fichier
    % csv il y a parfois 2 fois la même seconde répéété ce qui fait bugger
    % l'interpolation par la suite. On supprime ici les secondes en double
    datemat=datenum(str_date,'yyyy/mm/dd HH:MM:SS');
    I=find(diff(datemat)<=0);
    if length(I)>0
        P(I)=[];
        T(I)=[];
        S(I)=[];
        hrdata.F(I)=[];
        hrdata.O(I)=[];
        hrdata.L(I)=[];
        datemat(I)=[];
    end
    hrdata.date = datemat;      
    hrdata.P=P;
    hrdata.T=T;
    hrdata.S=S;
   

    % soucis dans les fichiers tres haute resolution la salinite a des
    % valeurs completements fausses a� la remontee ( svt egale a 0 ou negative)
    hrdata.S(hrdata.S<3)=NaN;
    
