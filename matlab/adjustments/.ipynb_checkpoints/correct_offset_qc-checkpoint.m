function correct_offset_qc(name_prof,smru_name,coeff,salinity_offsets)


%% correction de temperature/salinite base sur tests in situ + ORP
calib = ncload_struct(name_prof,'SCIENTIFIC_CALIB_COEFFICIENT');

Mqc.JULD = ncread(name_prof,'JULD');
I=1:length(Mqc.JULD);%find(strcmp(Mqc.platform_number,Mqc.list_descr{kk}));
if isempty(I), return; end

P1=0; P2=0; T1=0; T2=0; S1=0; S2=0; F1=0.6; F2=0;
list_var = {'T1','T2','S1','S2'};
for kk = 1:length(list_var),
    if ismember(smru_name,coeff.Properties.RowNames) & ...
            any(strcmp(list_var{kk},coeff(smru_name,:).Properties.VariableNames)) & ...
            coeff{smru_name,list_var{kk}} & ~isnan(coeff{smru_name,list_var{kk}})
        eval([list_var{kk} ' = coeff{smru_name,list_var{kk}};'])
    end
end

Mqc.PRES = ncread(name_prof,'PRES');
Mqc.PRES_ADJUSTED=Mqc.PRES;
Mqc.PRES_ADJUSTED(:,I)=Mqc.PRES(:,I) - P1 * 1e-3 * Mqc.PRES(:,I) - P2;
PRES2 = Mqc.PRES_ADJUSTED(:,I);
str_calib=sprintf('p1= %8.6e dbar/km, p2= %8.6e dbar',P1,P2);
str_calib=sprintf('%-256s',str_calib);
calib.SCIENTIFIC_CALIB_COEFFICIENT(:,1,1,I)=repmat(str_calib',1,length(I));

Mqc.TEMP = ncread(name_prof,'TEMP');
Mqc.TEMP_ADJUSTED=Mqc.TEMP;
Mqc.TEMP_ADJUSTED(:,I)=Mqc.TEMP(:,I) - T1 * 1e-3 * PRES2 - T2;
str_calib=sprintf('t1= %8.6e degC/km, t2= %8.6e degC',T1,T2);
str_calib=sprintf('%-256s',str_calib);
calib.SCIENTIFIC_CALIB_COEFFICIENT(:,2,1,I)=repmat(str_calib',1,length(I));

Mqc.PSAL = ncread(name_prof,'PSAL');
Mqc.PSAL_ADJUSTED=Mqc.PSAL;
Mqc.PSAL_ADJUSTED(:,I)=Mqc.PSAL(:,I) - S1 * 1e-3 * PRES2 - S2;
str_calib=sprintf('s1= %8.6e  psu/km, s2= %8.6e  psu',S1,S2);
str_calib=sprintf('%-256s',str_calib);
calib.SCIENTIFIC_CALIB_COEFFICIENT(:,3,1,I)=repmat(str_calib',1,length(I));

try, Mqc.CHLA = ncread(name_prof,'CHLA'); end
if isfield(Mqc,'CHLA')
    Mqc.CHLA_ADJUSTED = Mqc.CHLA;
    Mqc.CHLA_ADJUSTED(:,I) = F1 * Mqc.CHLA(:,I)+ F2;
    str_calib=sprintf('f1= 0.6, f2= 0.0');
    str_calib=sprintf('%-256s',str_calib);
    calib.SCIENTIFIC_CALIB_COEFFICIENT(:,4,1,I)=repmat(str_calib',1,length(I));
end

try, Mqc.DOXY = ncread(name_prof,'DOXY'); end
try, Mqc.LIGHT = ncread(name_prof,'LIGHT'); end

sal_cor = load_salinity_offset(smru_name,salinity_offsets,length(I));
if length(sal_cor)
    for ii=1:length(I)
        Mqc.PSAL_ADJUSTED(:,I(ii))=Mqc.PSAL_ADJUSTED(:,I(ii)) + sal_cor(ii);
        str_calib=sprintf('s1= %8.6e  psu/km, s2= %8.6e  psu', S1, S2 - sal_cor(ii));
        str_calib=sprintf('%-256s',str_calib);
        calib.SCIENTIFIC_CALIB_COEFFICIENT(:,3,1,I(ii)) = str_calib;
    end
end


%%

ncwrite(name_prof,'PRES_ADJUSTED',Mqc.PRES_ADJUSTED);
ncwrite(name_prof,'TEMP_ADJUSTED',Mqc.TEMP_ADJUSTED);
ncwrite(name_prof,'PSAL_ADJUSTED',Mqc.PSAL_ADJUSTED);

if isfield(Mqc,'CHLA')
    ncwrite(name_prof,'CHLA_ADJUSTED',Mqc.CHLA_ADJUSTED);
end

if isfield(Mqc,'DOXY')
    ncwrite(name_prof,'DOXY_ADJUSTED',Mqc.DOXY);
end

if isfield(Mqc,'LIGHT')
    ncwrite(name_prof,'LIGHT_ADJUSTED',Mqc.LIGHT);
end

ncwrite(name_prof,'SCIENTIFIC_CALIB_COEFFICIENT',calib.SCIENTIFIC_CALIB_COEFFICIENT);

