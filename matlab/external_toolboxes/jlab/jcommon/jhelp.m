function[varargout]=jhelp(str)
%JHELP  Opens linked JLAB help files in Matlab's internal web browser.
%
%   JHELP with no input arguments opens the JLAB Contents.m file.
%   JHELP MENU opens an index of all JLAB functions.
%   JHELP FILE opens the help information for JLAB function FILE.
%   
%   Usage: jhelp
%          jhelp menu
%          jhelp polysmooth
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2014--2015 J.M. Lilly --- type 'help jlab_license' for details

if nargin==0
    str='Contents';
end

%subdir=whichdir([str '.m']);
jlabdir=whichdir('jlab_license.m');
%subdir=subdir(length(jlabdir)+2:end);
%web([jlabdir '/doc/jlab/' subdir '/' str '.html'])
web([jlabdir '/doc/jlab/' str '.html'])

