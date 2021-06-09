%JLAB_ADDPATH  Adds JLAB subdirectories to your Matlab search path. 
%
%   This script adds JLAB subdirectories to your Matlab search path.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2015 J.M. Lilly --- type 'help jlab_license' for details

%First, find the path to JLAB
fullpath=which('jlab_addpath');
jlab_path=fullpath(1:end-15);

dirlist=dir(jlab_path);
for i=1:length(dirlist)
    if dirlist(i).isdir
        if ~strcmpi(dirlist(i).name(1),'.')
            addpath([jlab_path '/' dirlist(i).name])
        end
    end
end

clear i dirlist jlab_path fullpath