function ridgepack_globeatt(ncid,nc)

% ridgepack_globeatt - Adds global CF attributes to a netcdf file
%
% function ncglobal(ncid,nc)
%
% Adds CF-1.4 convention attributes to a netcdf file.  Inputs required are:
%
% INPUT:
%
% ncid    - ID number of the open netcdf file
% nc      - netcdf structure (defined below, but also see ridgepack_struct).
%
%
% OUTPUT:
%
% Appends to an open netcdf file in the working directory.
%
% Global attributes for the netcdf file appear in the nc data structure:
%
% nc.attributes.title - A succinct description of what is in the dataset.
%                       This is compulsory.
% nc.attributes.source - The method of production of the original data. 
%           If it is model-generated, source should name the model and its,
%           version as specifically as is useful. If it is an observational, 
%           source should characterize this (e.g., "surface observation" or 
%           "radiosonde"). This attribute should also include specifics of
%           of the instrumentation used to produce the data:
%           1) satellite/remotely sensed data: include summary 
%              of the satellite, instrument and processing algorithms.
%           2) GPS/positional data: include processing procedure and 
%              retrieval method.
%           3) directly sensed data: exact details on sensors used (model
%              numbers, sensor types)
%           4) model generated data: computer on which the data was 
%              generated.
%           5) conglomerate dataset (such as climatology or assimilated):
%              list (within reasonable bounds) all sources of data and 
%              methods/models used to compile the dataset.
%           This is compulsory.
% nc.attributes.references - Published or web-based references that describe 
%           the data or methods used to produce it. This need to be included 
%           if references are not applicable.  Exclusion of this will
%           will remove the entry from the netcdf file. This is optional.
% nc.attributes.comment - Miscellaneous information about the data or
%           methods use to produce it.  This is optional.
%
% If attribute text is long, it is advisable to place newline markers 
% at appropriate placed within the text (in netcdf files these are "\n").
%
% This function also requires two environment variables to be supplied
% in the UNIX/Linux shell for the global attributes when it is run: 
%
% USER_ORG - the institution or organization in which the file is created.
% USER_REAL - the real name of the user.  
%
% These are supplied by environment variable rather than direct input
% to improve portability of the function, and to remove replication input.
% They may be set in bash using: export USER_REAL='Joe Blogs', or in csh 
% using: setenv USER_REAL 'Joe Blogs'. These commands can be added to 
% .bashrc or .cshrc, files respectively. They may also be set in the
% startup.m file for matlab using, e.g.:
%
% setenv('USER_ORG','International Arctic Research Center');
% setenv('USER_REAL','Andrew Roberts');
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if isfield(nc,'attributes') & ~isstruct(nc.attributes) ;
        error('attributes input data is not a structure');
end

% set up basic global attributes
if (isfield(nc,'attributes') & ~isfield(nc.attributes,{'title','TITLE','Title'})) | ~isfield(nc,'attributes')
	nc.attributes.title=input('Enter the netcdf file title: ','s');
end

if ~isfield(nc.attributes,{'source','SOURCE','Source'})
	nc.attributes.source=[getenv('HOST'),':',pwd];
end

if (debug & ~isfield(nc.attributes,{'references','REFERENCES','References'})); 
 	disp(['No references attributed to netcdf file'])
end

if (debug & ~isfield(nc.attributes,{'comment','COMMMENT','Comment'}));
 	disp(['No comments attributed to netcdf file'])
end

if isfield(nc.attributes,{'conventions'});
	nc.attributes=rmfield(nc.attributes,'conventions');
end
if isfield(nc.attributes,{'CONVENTIONS'});
	nc.attributes=rmfield(nc.attributes,'CONVENTIONS');
end
nc.attributes.Conventions='CF-1.4';


% Obtain file history and institution attribute information from the shell
if isunix | ismac;
	user_org=getenv('USER_ORG');
	if isempty(user_org);
		user_org=input('Enter your affiliation: ','s'); 
		setenv('USER_ORG',user_org);
        end
	user_real=getenv('USER_REAL');
	if isempty(user_real);
		user_org=input('Enter your name: ','s'); 
		setenv('USER_REAL',user_real);
        end
else
	S=license('inuse','MATLAB');
	user_real=char(S.user);
	user_org=[];
end

if  ~isfield(nc.attributes,'history')
  nc.attributes.history=[datestr(now),': File created by ',user_real];
end

if  ~isfield(nc.attributes,'institution') & isempty(user_org);
  disp(['No institution attributed to netcdf file'])
elseif  ~isfield(nc.attributes,'institution');
  nc.attributes.institution=user_org;
end


% Assign global attributes
gatt=fieldnames(nc.attributes);
gid=netcdf.getConstant('NC_GLOBAL');
netcdf.reDef(ncid);
for i=1:length(gatt)
  attname=char(gatt{i});
  if ischar(nc.attributes.(attname)) & ~strcmp(nc.attributes.(attname),'')
    netcdf.putAtt(ncid,gid,attname,char(nc.attributes.(attname)));
  else
    disp([attname,' is not a string and will not be written to global attributes.']);
  end
end
netcdf.endDef(ncid);

if debug; disp(['...Leaving ',mfilename]); end

