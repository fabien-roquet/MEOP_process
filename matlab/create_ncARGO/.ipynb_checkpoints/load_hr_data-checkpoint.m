function hrdata = load_hr_data(name_hr_file,continuous)


    hrdata = [];
    
    %% chargement du fichier ctd haute resolution
    fid = fopen(name_hr_file);
    tline = fgetl(fid);
    fclose(fid);
    
    Fluo = [strfind(tline,'Fluo') strfind(tline,'FLUORO')];
    Oxy = strfind(tline,'Oxy');
    Light = strfind(tline,'PPFD');
    hrdata.isfluo=length(Fluo); hrdata.isoxy=length(Oxy); hrdata.islight=length(Light);
    hrdata.continuous = continuous;
    
    try
        if ~hrdata.isfluo & ~hrdata.isoxy & ~hrdata.islight
            [str_date,P,T,S] = ...
                textread(name_hr_file,'%s%f%f%*f%f%*[^\n]',...
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
    
    hrdata.date = datenum(str_date,'yyyy/mm/dd HH:MM:SS');      
    hrdata.P=P;
    hrdata.T=T;
    hrdata.S=S;

    % soucis dans les fichiers tres haute resolution la salinite a des
    % valeurs completements fausses a� la remontee ( svt egale a 0 ou negative)
    hrdata.S(hrdata.S<3)=NaN;
    
