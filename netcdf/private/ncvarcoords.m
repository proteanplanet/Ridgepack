function [coords]=ncvarcoords(nc,var)

% NCVARCOORDS - Lists the coordinates in a cell array of a given variable 
%
% function [coords]=ncvarcoords(nc,var)
%
% This function takes apart the space-separated coordinates set in the netcdf
% structure nc and turns them into a cell array.  Coordinate variables are 
% recognized with the coordinate attribute in the CF convention such as:
%
% nc.temperature.coordinate='latitude longitude' 
%
% where each coordinate variable is separated by a space and the attribute
% is a simple character array.  
%
% Input:
% nc  - netcdf structure (see ncstruct)
% var - given variable in a netcdf structure
% 
% Output:
% coords - coordinates in a cell array for the given variable
%
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
%
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% run basic checks
if ~isstruct(nc)
	error('Input to ncvarcoords is not a structure')
elseif ~isfield(nc,var)
	error([var,' is not a variable in nc'])
end

coords={};

if isfield(nc.(var),'coordinates')

 if iscell(nc.(var).coordinates)
  error([var,' coordinates is set as a cell array']);
 end

 j=0;

 [token,remain]=strtok(nc.(var).coordinates);

 if ~any(strcmp(token,coords))
	 j=j+1;
	 coords{j}=token;

	 % check that token is a variable
	 if ~isfield(nc,token) & debug
 	   disp(['WARNING: The coordinate ',token,' listed as a coordinate for ',...
	          var,' is not a variable']);
  	 end

 end

 while ~isempty(remain)

	 [token,remain]=strtok(remain);
			 
	 if ~any(strcmp(token,coords))
		 j=j+1;
		 coords{j}=token;

	 	 % check that token is a variable
	 	 if ~isfield(nc,token) & debug
		   disp(['WARNING: The coordinate ',token,' listed as a coordinate for ',...
		           var,' is not a variable']);
  	 	 end

	 end

 end

end

if debug; disp(['Coordinates set for ',var]); end

if debug; disp(['...Leaving ',mfilename]); end

