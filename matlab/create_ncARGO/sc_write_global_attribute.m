%% script updating metadata (global attributes)

ncwriteatt(name_prof,'/','pi_name',info_deployment.PI);
ncwriteatt(name_prof,'/','data_type','Marine mammals time-series data');
ncwriteatt(name_prof,'/','format_version','1.1');
ncwriteatt(name_prof,'/','date_update',datestr(now(),'yyyy-mm-ddTHH:MM:SSZ'));
ncwriteatt(name_prof,'/','version_database',conf.version);
ncwriteatt(name_prof,'/','PI',info_deployment.PI);
ncwriteatt(name_prof,'/','reference_file_name',sprintf('%s_prof.nc',smru_name));
if info_deployment.public,
    ncwriteatt(name_prof,'/','is_the_data_public','yes: data can be used freely providing that data owner is properly cited (see meop.net for citing information)');
else
    ncwriteatt(name_prof,'/','is_the_data_public','no: data cannot be used without the prior consent of the data owner');
end
ncwriteatt(name_prof,'/','nation',info_deployment.NATION);
ncwriteatt(name_prof,'/','deployment_code',info_deployment.EXP);
ncwriteatt(name_prof,'/','source','Marine mammal observation');
ncwriteatt(name_prof,'/','data_mode','D');
ncwriteatt(name_prof,'/','references','http://www.meop.net');
ncwriteatt(name_prof,'/','reference_doi',' ');
ncwriteatt(name_prof,'/','Conventions','CF-1.6 Sea-mammals-1.1');
ncwriteatt(name_prof,'/','Netcdf_version','3.6');
ncwriteatt(name_prof,'/','naming_authority','MEOP consortium (Marine Mammals Exploring the Oceans Pole to Pole)');
ncwriteatt(name_prof,'/','cdm_data_type','Station');
ncwriteatt(name_prof,'/','geospatial_vertical_min',0);
ncwriteatt(name_prof,'/','geospatial_vertical_max',2000);
ncwriteatt(name_prof,'/','data_assembly_center','MEOP/Fabien Roquet (MISU)');
ncwriteatt(name_prof,'/','distribution_statement','Follow MEOP data policy standards, cf. http://www.meop.net/the-dataset/data-access.html. Data available free of charge. User assumes all risk for use of data. User must display citation in any publication or product using data. User must contact PI prior to any commercial use of data');
ncwriteatt(name_prof,'/','citation','The marine mammal data were collected and made freely available by the International MEOP Consortium and the national programs that contribute to it (http://www.meop.net).');
ncwriteatt(name_prof,'/','thermal_lag_adjustment','no');

[smru_prefix,Nsplit] = Nsplit_from_smru_name(smru_name);
K=find(strcmp(conf.list_smru_platform_code,smru_prefix));
platform = conf.platform_json(K);
try
    ncwriteatt(name_prof,'/','platform_code',platform.platform_code);
    if length(platform.wmo_platform_code)>0
        ncwriteatt(name_prof,'/',...
            'wmo_platform_code',sprintf('Q99%05d',str2num(platform.wmo_platform_code)));
    else
        ncwriteatt(name_prof,'/','wmo_platform_code',' ');
    end
    ncwriteatt(name_prof,'/','smru_platform_code',smru_name);
    ncwriteatt(name_prof,'/','deployment_code',platform.deployment_code);
    ncwriteatt(name_prof,'/','species',platform.species);
    ncwriteatt(name_prof,'/','time_coverage_start',platform.time_coverage_start);
    if isempty(platform.time_coverage_end)
        if length(strfind(name_prof,'traj'))>0
            juld=ncread(name_prof,'TIME');
        else
            juld = ncread(name_prof,'JULD');
        end
        juld_ref = datenum('19500000','yyyymmdd');
        if isempty(juld),
            platform.time_coverage_end = platform.time_coverage_start;
        else
            platform.time_coverage_end = datestr(max(juld)+juld_ref,'yyyy-mm-ddT00:00:00Z');
        end
    end
    ncwriteatt(name_prof,'/','time_coverage_end',platform.time_coverage_end);
    location = platform.location; location(strfind(location,','))='.';
    ncwriteatt(name_prof,'/','location',location);
    if isfield(platform,'loc_algorithm')
        ncwriteatt(name_prof,'/','loc_algorithm',platform.loc_algorithm);
    end
    ncwriteatt(name_prof,'/','firmware_version',platform.firmware_version);
    ncwriteatt(name_prof,'/','firmware_parameters',platform.firmware_parameters);
    ncwriteatt(name_prof,'/','instr_id',platform.instr_id);
    ncwriteatt(name_prof,'/','ptt',platform.ptt);
catch
    disp(['problem encountered while writing metadata in the ncARGO file :' smru_name])
end

%% manual update of metadata
[smru_prefix,Nsplit] = Nsplit_from_smru_name(smru_name);
K2=find(ismember(conf.table_meta.smru_platform_code,smru_prefix));
if length(K2)==1,
    var = conf.table_meta.Properties.VariableNames;
    for kk=1:length(var)
        if iscell(conf.table_meta{K2,kk}) && length(conf.table_meta{K2,kk}{:}),
            ncwriteatt(name_prof,'/',var{kk},conf.table_meta{K2,kk}{:});
        end
    end
elseif length(K)>1,
    disp(['problem: multiple entry for ' smru_name ' in table_meta.csv']);
end




