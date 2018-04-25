function [nc]=ridgepack_add(nc,name,data,long_name,dimension,units,type,title)

% ridgepack_add - Add data and descriptive metadata to a netCDF structure
%
% function [nc]=ridgepack_add(nc,name,data,long_name,dimension,units,type,title)
%
% This function adds a variable to a netCDF structure.  To learn
% more about netCDF structures, see function ridgepack_struct.  The nc structure
% variable nc may be omitted from the input list if the structure
% does not exist, so that the function can also be run as:
%
% function [nc]=ridgepack_add(name,data,long_name,dimension,units,type,title)
%
% INPUT:
%
% nc        - netCDF structure (can be ommitted for starting a new structure)
% name      - name of variable to be added
% data      - data to be added
% long_name - long name to be used in netCDF file
% dimension - cell array of relevant dimensions
% units     - units of data (optional or enter '' to leave blank)
% type      - netCDF type (e.g.'NC_CHAR','NC_DOUBLE')
% title     - optional character arguement giving the title of the nc structure
%
%
% OUTPUT:
%
% nc        - netCDF structure with added data
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if ischar(nc);
 if nargin>6; title=type; end
 if nargin>5; type=units; end
 if nargin>4; units=dimension; end
 if nargin<4 ; error('Not enough arguments'); end
 dimension=long_name;
 long_name=data;
 data=name;
 name=nc;
 clear nc;
 if nargin > 6 & ischar(title);
  nc.attributes.title=title;
 else
  nc.attributes.title=input('Enter the title of the nc structure: ','s');
 end
elseif ~isstruct(nc)
 error('nc is incorrectly specified');
elseif nargin<5 ;
 error('Not enough arguments')
end

nc.(name).data=data;

if ischar(long_name);
	nc.(name).long_name=long_name;
else
	error('long_name must be a character variable');
end


if iscell(dimension);
	nc.(name).dimension=dimension;
else
	error('dimension must be a cell array')
end


if nargin > 5 & ischar(units) & ~strcmp(units,'');
	nc.(name).units=units;
elseif nargin > 5 & ~strcmp(units,'');
	error('units must be a character variable');
end

if nargin > 6 & ischar(type);
	nc.(name).type=type;
else
	nc.(name).type='NC_DOUBLE';
end

if debug; disp(['...Leaving ',mfilename]); end
