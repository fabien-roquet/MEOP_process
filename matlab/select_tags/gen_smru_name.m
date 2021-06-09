function smru_name = gen_name_prof(smru_prefix,Nsplit)


    if length(Nsplit) & isnumeric(Nsplit)
        Nsplit=['-N' num2str(Nsplit)];
    elseif length(Nsplit) & ~strcmp(Nsplit(1),'-')
        Nsplit=['-' Nsplit];
    end
    
    smru_name = sprintf('%s%s',smru_prefix,Nsplit);
    
