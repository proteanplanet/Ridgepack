% This is a sample MATLAB startup file to accompany Ridgepack 
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)

% Check this is a unix system
if not(isunix); 
 error('This matlab setup designed for Unix/Linux/MacOSX based systems'); 
end

% Set paths
path(genpath([getenv('HOME'),'/local/matlab']),path)
path(genpath([getenv('HOME'),'/Ridgepack']),path)

% Defaults fonts and interpreter
set(0,'DefaultTextInterpreter','Latex')
set(0,'defaultaxesfontname','Serif')
set(0,'defaulttextfontname','Serif')

% Set defaults and environment variable for Ridgepack functions
%  organization and name that appears in netcdf files written
setenv('USER_ORG','my organization'); 
setenv('USER_REAL','My Name');
%  location of data files (ETOPO and WRF files for RASM)
setenv('ETOPOFILE',[getenv('HOME'),'/data/TOPOBATHYMETRY/ETOPO2/ETOPO2.raw.bin'])
setenv('WRFGRIDFILE',[getenv('HOME'),'/data/MODEL/RACM'])

% Add further unix paths, including Fink if using Mac OS X
unixdirs={'/usr/local/bin','/sw/bin','/sw/sbin'};
for i=1:length(unixdirs)
 dirs=char(unixdirs{i});
 A=exist(dirs);
 B=findstr(getenv('PATH'),dirs);
 if A==7 & isempty(B); setenv('PATH',[getenv('PATH'),':',dirs]); end
end
clear A B dirs unixdirs i

% basic display
format LONGG;
format COMPACT;

% Set default directory
if exist([getenv('HOME'),'/work'],'dir')>0
 cd([getenv('HOME'),'/work']);
else
 disp('Work directory is missing from home directory')
end

