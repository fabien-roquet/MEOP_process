function [smru_prefix,Nsplit] = Nsplit_from_smru_name(smru_name)

    smru_prefix = smru_name;
    Nsplit = '';
    if length(smru_name)>3 & strcmp(smru_name(end-2:end-1),'-N'), % split tag
        smru_prefix = smru_name(1:end-3);
        Nsplit = smru_name(end-1:end);        
    end    
    