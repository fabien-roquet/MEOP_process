function [smru_prefix,Nsplit,suffix,namedir] = smru_name_from_name_prof(name_prof)

    if length(name_prof)<12,
        error(['Wrong format of name_prof variable: ' name_prof])
    end
    
    namedir='';
    I=strfind(name_prof,'/');
    if length(I),
        namedir = name_prof(1:I(end));
    else
        I=0;
    end
    
    smru_name = name_prof(I(end)+1:end-12);
    [smru_prefix,Nsplit] = Nsplit_from_smru_name(smru_name);
    suffix = name_prof(end-10:end-8);
    
