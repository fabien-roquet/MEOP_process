function [offset]=load_salinity_offset(smru_name,salinity_offsets,N_profiles)

use_salinity_offset = any(strcmp(smru_name,salinity_offsets.Properties.RowNames));
if use_salinity_offset
    
    coeffs = salinity_offsets(smru_name,:);
    index = []; offset = [];
    for kk = 1:4,
        index(kk) = coeffs{smru_name,sprintf('index_%d',kk)};
        offset(kk) = coeffs{smru_name,sprintf('offset_%d',kk)};
    end
    
    J = find(index==0);
    index(J(1)) = N_profiles;
    
    if length(J)>1,
        index (J(2:end)) = [];
        offset(J(2:end)) = [];
    end
    
    offset= interp1(index,offset,1:N_profiles);
    
else
    
    offset = [];
    
end
