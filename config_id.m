function config_id = config_id()
% determine an identifier of the current machine

    if ispc, name = getenv('COMPUTERNAME');
    else
        [~, name] = system('hostname');
        name = [getenv('USER') '_' deblank(name) '_linux'];
    end
    name = strtrim(lower(name));
    config_id = strrep(name,'.','_');    
    config_id = strrep(config_id,'-','_');    
