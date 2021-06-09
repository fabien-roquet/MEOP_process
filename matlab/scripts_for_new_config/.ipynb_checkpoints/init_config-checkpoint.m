function conf = init_config(new_processdir)
% set up matlab environment and initialize it if first time it is used
% The first time you want to install a config:
%  1) copy files from the folder scripts_for_new_config in the folder
%       where you want to process data,
%  2) set folder paths in the json file configs.json. 
%       Run the command config_id to obtain the config name you must
%       use for your computer. Enter absolute path for:
%           "processdir":"xxx", % path where to run processing.
%           "datadir":"xxx", % raw data
%           "matlabdir":"?/GitHub/meop_qc/matlab_toolbox/",
%           "inputdir":"?/", % SMRU input data
%           "refdir":"?/", % folder with reference climatologies
%           "public":"?/", % folder to write public data
%           "pdflatex":"pdflatex" %command to run the pdflatex compiler
% IMPORTANT: To install the new config, go to the folder scripts_for_new_config,
%  set configs.json and provide in particular the path to intended processdir folder.
%  Once run, move to the processdir folder and from there run the process scripts.

%% initialization of the mirounga processing system
d = jsondecode(fileread('configs.json'));
if isempty(d.config),
    d.config = config_id();    
end
try
    conf = getfield(d.configs,d.config);
catch
    error([d.config ' is not a supported config. modify configs.json']);
end
conf.version     = d.version.CTDnew;
conf.version_old = d.version.CTDold;
conf.version_SMS = d.version.SMSnew;

% set up matlab path
p=genpath(conf.matlabdir);
Ibeg=[1 strfind(p,':')+1]; Ibeg(end)=[];
Iend=[strfind(p,':')-1];
if length(Ibeg),
    p2=[];
    for ii=1:length(Ibeg),
        name=p(Ibeg(ii):Iend(ii));
        if ~length(strfind(name,'/.AppleDouble'))
            p2=[p2 ':' name ];
        end
    end
    p2(1)=[]; p2(end+1)=':';
end
addpath(p2)
rmpath([conf.matlabdir 'scripts_for_new_config'])

conf.rawdir          = [conf.datadir 'raw_smru_data_odv/'];
conf.rawdir_hr       = [conf.datadir 'raw_smru_hr_data/'];
conf.json            = [conf.datadir 'config_files/'];
conf.crawl.locdir    = [conf.datadir 'crawl_locations/'];
conf.cls.locdir      = [conf.datadir 'smooth_cls_locations/'];

conf.datadir         = [conf.processdir 'final_dataset_prof/'];
conf.mapsdir         = [conf.processdir 'maps/'];
conf.texdir          = [conf.processdir 'doc_latex/'];
conf.plotdir         = [conf.processdir 'plots/'];
conf.calibplotdir    = [conf.processdir 'calibration_plots/'];
conf.csv_config      = [conf.processdir 'csv_config_files/'];

conf.woddir          = [conf.refdir 'WOD_data/WOD_matlab_nc/'];
conf.coradir         = [conf.refdir 'CORA_data/CORA_ncfiles/'];
conf.meopdir         = [conf.refdir 'MEOP_last_stable_version/'];

conf.temporary       = [conf.processdir 'temporary/'];
conf.temporary_tex   = [conf.processdir 'temporary/tex/'];
conf.temporary_fcell = [conf.processdir 'temporary/fcell/'];

% create process_dir if not already existing
if ~exist(conf.processdir,'dir'),
    try,
        [s,mess,messid] = mkdir(conf.processdir);
        [s,mess,messid] = copyfile([conf.matlabdir 'csv_config_files_ref/'],conf.processdir);
        [s,mess,messid] = copyfile([conf.matlabdir 'init_config.m'],conf.processdir);
        [s,mess,messid] = mkdir(conf.datadir);
        [s,mess,messid] = mkdir(conf.mapsdir);
        [s,mess,messid] = mkdir([conf.mapsdir 'deployments']);
        [s,mess,messid] = mkdir([conf.mapsdir 'groups']);
        [s,mess,messid] = mkdir([conf.mapsdir 'global']);
        [s,mess,messid] = mkdir(conf.texdir);
        [s,mess,messid] = mkdir(conf.plotdir);
        [s,mess,messid] = mkdir(conf.calibplotdir);
        [s,mess,messid] = mkdir(conf.temporary);
        [s,mess,messid] = mkdir(conf.temporary_tex);
        [s,mess,messid] = mkdir(conf.temporary_fcell);
        disp(['Processing done in: ' conf.processdir])
        cd(conf.processdir);
    catch
        delete(conf.processdir);
    end
else
    cd(conf.processdir);
end



