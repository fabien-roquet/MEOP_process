function[varargout]=jlab_allhelp(varargin)
%JLAB_ALLHELP  Displays the help comments for all JLAB modules.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details
 
jlab_path=findpath('jlab');
dirlist=dir(jlab_path);
for i=1:length(dirlist)
    if dirlist(i).isdir
        if ~strcmpi(dirlist(i).name(1),'.')
            if ~strcmpi(dirlist(i).name,'doc')&&~strcmpi(dirlist(i).name,'html')&&~strcmpi(dirlist(i).name,'figures')
                eval(['help ' dirlist(i).name])
            end
        end
    end
end