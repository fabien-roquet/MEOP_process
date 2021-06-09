function table_coeff = table_coeff_read(conf,EXPorTAG)
% table_coeff = table_coeff_read(conf,EXPorTAG)
%   conf is a conf struct obtained from init_mirounga. Can be empty.
%   EXPorTAG is a string. All entries with names starting with EXPorTAG will be returned.

if isempty(conf),
    conf = init_mirounga;
end

if ~exist('EXPorTAG','var') % all tags from EXP deployment
    table_coeff = conf.table_coeff;
else
    table_coeff = conf.table_coeff(startsWith(conf.table_coeff.Properties.RowNames,EXPorTAG),:);
end

