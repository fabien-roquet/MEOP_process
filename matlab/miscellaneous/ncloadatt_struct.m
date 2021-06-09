function data = ncloadatt_struct(fileIn, varargin);

% ncload -- Load NetCDF variables.
%  ncload('fileIn', 'var1', 'var2', ...) loads the
%   given variables of 'fileIn' into the Matlab
%   workspace of the "caller" of this routine.  If no names
%   are given, all variables are loaded.

data=struct;

if ~exist(fileIn,'file'), 
    disp(sprintf('error: %s does not exist',fileIn))
    return,
end

if isempty(varargin); 
    finfo = ncinfo(fileIn);
    atts={finfo.Attributes.Name};
    varargin = atts; 
end;

for ii = 1:length(varargin)
    try
        aa = ncreadatt(fileIn,'/',varargin{ii});
        if isnumeric(aa), aa=double(aa); end
        data = setfield( data, varargin{ii}, aa );
    catch
        disp(sprintf('error: %s is not an attribute of %s',varargin{ii},fileIn));
    end;
end


