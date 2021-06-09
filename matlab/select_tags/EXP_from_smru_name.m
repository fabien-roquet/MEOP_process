function EXP = EXP_from_smru_name(smru_name)

    I = strfind(smru_name,'-');
    EXP = smru_name(1:I(1)-1);
    
    