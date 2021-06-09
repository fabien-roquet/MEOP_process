function table_coeff_write(conf,new_table_coeff)
% table_coeff_write(conf,table_coeff)
%   conf is a conf struct obtained from init_mirounga. Can be empty.
%   table_coeff is a table with the coefficients that will be written in the table_coeff.csv file.

if isempty(conf),
    conf = init_mirounga;
end

table_coeff = conf.table_coeff;

for kk=1:length(new_table_coeff.Properties.RowNames),
    [LIA,LOCB] = ismember(new_table_coeff.Properties.RowNames{kk},table_coeff.Properties.RowNames)
    if LIA,
        table_coeff(LOCB,:)=new_table_coeff(kk,:);
    end
end

name_file=[conf.processdir 'table_coeff.csv'];
if ~exist(name_file,'file')
    error(['WARNING: the file ' name_file ' was not found!'])
end
writetable(table_coeff,name_file,'WriteRowNames',1,'Delimiter',',');

