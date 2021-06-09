% script to process new deployments
% run deployments listed according to select_case, except if already done (see list 
%    in list_file_done)
% In addition, deployments in list_forced will be done regardless of previous criteria
init_config;

select_case = 1;
list_file_done = 'list_deployments_DONE.txt';
list_forced={};%'ct144','ct148','ct149','ft24','ft22','ct152'};

switch select_case,
    case 1
        % process all deployments to be processed
        EXPs = tags_processed;
        lEXP = {EXPs.Properties.RowNames{:}};
    case 2
        % process selected deployments to be processed
        EXPs = tags_processed([],list_forced);
        lEXP = {EXPs.Properties.RowNames{:}};
    case 3
        % select all deployments from a NATION to be processed
        EXPs = tags_processed([],'FRANCE');
        lEXP = {EXPs.Properties.RowNames{:}};
    case 4
        % process all deployments from a NATION
        list_file = 'list_deployments.txt';
        lEXP = readtable(list_file,'Delimiter',',','ReadVariableNames',0);
        lEXP = {lEXP{:,1}{:}};
    otherwise
        lEXP = {};
end

% add deployments from list_forced if not already present
for kk=1:length(list_forced),
    if ~any(ismember(lEXP,list_forced{kk})),
        lEXP = {list_forced{kk},lEXP{:}};
    end
end

% read list of files already done
lEXP_done = {};
if exist(list_file_done,'file')
    lEXP_done = readtable(list_file_done,'Delimiter',',','ReadVariableNames',0);
    lEXP_done = {lEXP_done{:,1}{:}};
end

% main loop
error=0;
for kk = 1:length(lEXP),
    EXP = lEXP{kk};
    if ismember(EXP,lEXP_done) & ~ismember(EXP, list_forced),
        disp([EXP ' already done: skipped']);
        continue
    end
    disp(' ')
    disp(EXP)
    try
        process_single_deployment(EXP);
        if ~ismember(EXP,lEXP_done)
            fid2 = fopen(list_file_done,'a');
            data = fprintf(fid2, '%s,\n', EXP);
            fclose(fid2);
        end
    catch exception
        msgText = getReport(exception,'extended','hyperlinks','off');
        fprintf(1,'%s',msgText);
        error = 1;
    end
end

if ~error
    % plot_global_dataset(conf);
    % plot_global_dataset_SMS(conf);    
    % sc_load_data_to_public_folder;
    % sc_load_data_to_public_folder_SMS;
end


