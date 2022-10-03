  
function fusion_profilTS_profilFL(EXP,path)

namefile = [path  '/' EXP '_ODV.txt'];
namefile_FL = [path  '/' EXP '_FL_ODV.txt'];
namefile_CTD = [path  '/' EXP '_CTD_ODV.txt'];

copyfile(namefile,namefile_CTD);

fid = fopen(namefile_CTD);
tline = fgetl(fid);
tline = fgetl(fid);
fclose(fid);

[tag,numProf,type,date,lon,lat,bot,P,T,S,C] = ...
        textread(namefile_CTD,'%s%d%s%s%f%f%d%f%f%f%f',...
        'delimiter',';','headerlines',2);
   
  

fid = fopen(namefile_FL);
allText = textscan(fid,'%s','delimiter','\n');
numberOfLines = length(allText{1});
fclose(fid);

if numberOfLines < 3,
   return
end

fid = fopen(namefile_FL);
tline = fgetl(fid);
tline = fgetl(fid);
fclose(fid);

[tagfl,numProffl,datefl,lonfl,latfl,Pfl,F,L,O] = ...
        textread(namefile_FL,'%s%d%*s%s%f%f%*d%f%f%f%f%*f',...
        'delimiter',';','headerlines',2);
        
        
%%
Obis=[];
Fbis=[];
Lbis=[];
tagbis=[];
numProfbis=[];
datebis=[];
lonbis=[];
latbis=[];
Pbis=[];
Tbis=[];
Sbis=[];
typebis=[];
botbis=[];
A=[date,tag];
[profil,ia,ic]=unique(A(:,1));

for ii=1:length(A)
    if length(A{ii,1})>0
        Ideb=ii;
        Ifin=find(numProf(Ideb+1:end)~=0);
        if Ifin
            Ifin=Ideb+Ifin(1)-1;
        else
            Ifin=length(numProf);
        end
        start=find(strcmp(datefl,A(ii,1)) & strcmp(A(ii,2),tagfl));
        
        if start
            if length(start)>1
                Jdeb=start(1);
                Jdebbis=start(2);
                Jfin=find(numProffl(Jdeb+1:end)~=0);
                Jfinbis=find(numProffl(Jdebbis+1:end)~=0);
                if Jfin
                    Jfin=Jdeb+Jfin(1)-1;
                else
                    Jfin=length(numProffl);
                end
                if Jfinbis
                    Jfinbis=Jdebbis+Jfinbis(1)-1;
                else
                    Jfinbis=length(numProffl);
                end
                Prof=unique([Pfl(Jdeb:Jfin);P(Ideb:Ifin)]);
                if length(find(F(Jdeb:Jfin)~=999))>0
                    Fprof=interp1(Pfl(Jdeb:Jfin),F(Jdeb:Jfin),Prof);
                else
                    Fprof=interp1(Pfl(Jdebbis:Jfinbis),F(Jdebbis:Jfinbis),Prof);
                end
                if length(find(L(Jdeb:Jfin)~=999))>0
                    Lprof=interp1(Pfl(Jdeb:Jfin),L(Jdeb:Jfin),Prof);
                else
                    Lprof=interp1(Pfl(Jdebbis:Jfinbis),L(Jdebbis:Jfinbis),Prof);
                end
                Oprof=interp1(Pfl(Jdeb:Jfin),O(Jdeb:Jfin),Prof);
            else
                Jdeb=start;
                Jfin=find(numProffl(Jdeb+1:end)~=0);
                if Jfin
                    Jfin=Jdeb+Jfin(1)-1;
                else
                    Jfin=length(numProffl);
                end
                Prof=unique([Pfl(Jdeb:Jfin);P(Ideb:Ifin)]);
                
                Fprof=interp1(Pfl(Jdeb:Jfin),F(Jdeb:Jfin),Prof);
                
                Lprof=interp1(Pfl(Jdeb:Jfin),L(Jdeb:Jfin),Prof);
                
                Oprof=interp1(Pfl(Jdeb:Jfin),O(Jdeb:Jfin),Prof);
            end
                
                emptyCell = cell(1,length(Prof)-1)';
                tagbis=[tagbis;tag(Ideb);emptyCell];
                typebis=[typebis;type(Ideb);emptyCell];
                numProfbis=[numProfbis;numProf(Ideb);numProf(1:length(Prof)-1)*0];
                datebis=[datebis;date(Ideb);emptyCell];
                lonbis=[lonbis;lon(Ideb);numProf(1:length(Prof)-1)*0];
                latbis=[latbis;lat(Ideb);numProf(1:length(Prof)-1)*0];
                botbis=[botbis;bot(Ideb);numProf(1:length(Prof)-1)*0];
                Pbis=[Pbis;Prof];
                Fbis=[Fbis;Fprof];
                Lbis=[Lbis;Lprof];
                Obis=[Obis;Oprof];
                Tbis=[Tbis;interp1(P(Ideb:Ifin),T(Ideb:Ifin),Prof)];
                Sbis=[Sbis;interp1(P(Ideb:Ifin),S(Ideb:Ifin),Prof)];
            
        else
            Obis=[Obis;T(Ideb:Ifin)*0+999];
            Fbis=[Fbis;T(Ideb:Ifin)*0+999];
            Lbis=[Lbis;T(Ideb:Ifin)*0+999];
            tagbis=[tagbis;tag(Ideb:Ifin)];
            typebis=[typebis;type(Ideb:Ifin)];
            numProfbis=[numProfbis;numProf(Ideb:Ifin)];
            datebis=[datebis;date(Ideb:Ifin)];
            lonbis=[lonbis;lon(Ideb:Ifin)];
            latbis=[latbis;lat(Ideb:Ifin)];
            botbis=[botbis;bot(Ideb:Ifin)];
            Pbis=[Pbis;P(Ideb:Ifin)];
            Tbis=[Tbis;T(Ideb:Ifin)];
            Sbis=[Sbis;S(Ideb:Ifin)];
          
        end
    end
    
end
I=find(isnan(Fbis));
Fbis(I)=999;
I=find(isnan(Lbis));
Lbis(I)=999;
I=find(isnan(Obis));
Obis(I)=999;
I=find(isnan(Tbis));
Tbis(I)=999;
I=find(isnan(Sbis));
Sbis(I)=999;

fid=fopen(namefile,'w');
fprintf(fid,'// Generic ODV file \n'); 
fprintf(fid,'Cruise;Station;Type;yyyy-mm-dd hh:mm;Longitude [degrees_east];Latitude [degrees_north];Bot. Depth [m];Pressure [dbar];Temperature [C];Salinity [PSU];Fluorescence [mg/l];Light [ln(PPFD)];Oxygen \n');

for ii=1:length(tagbis)
    
    fprintf(fid,'%s;%d;%s;%s;%f;%f;%d;%f;%f;%f;%f;%f;%f \n',tagbis{ii},numProfbis(ii),typebis{ii},datestr(datebis{ii},'yyyy-mm-dd HH:MM'),lonbis(ii),latbis(ii),botbis(ii),Pbis(ii),Tbis(ii),Sbis(ii),Fbis(ii),Lbis(ii),Obis(ii));
end
fclose(fid);
