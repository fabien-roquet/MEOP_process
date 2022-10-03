function import_ODV_data(conf,deployment_code)
% function that looks for zip data in the smru folder (conf.inputdir)
% copy zip files EXP_ODV.zip and EXP_FL_ODV.zip (if exists) locally, unzip and merge
% add new version of EXP_ODV.txt file into conf.rawdir

unzip([conf.inputdir deployment_code '/' deployment_code '_ODV.zip'],[conf.inputdir deployment_code]);

if length(dir([conf.inputdir deployment_code '\*FL*']))>0
   fusion_profilTS_profilFL(deployment_code,[conf.inputdir deployment_code])
end

copyfile([conf.inputdir deployment_code '\' deployment_code '_ODV.txt'],[conf.rawdir '/'])

disp(sprintf('%s: import in datadir',deployment_code));
rmdir([conf.inputdir deployment_code],'s')






