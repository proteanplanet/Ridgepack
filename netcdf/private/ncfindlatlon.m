function [lonvar,latvar]=ncfindlatlon(nc,coor);

% NCFINDLATLON - Find latitude and longitude in subset of variable in nc structure
%
% function [lonvar,latvar]=ncfindlatlon(nc,coor);
%
% This function searches the cell array coor for latitude- and longitude-
% like variables in the nc structure, where coor is a subset of variable
% names in the netcdf structure nc.
%
% Input:
% nc   - netcdf structure (see ncstruct for details)
% coor - cell array containing a subset of variable names in nc
%
% Output:
% lonvar - string of the longitude coordinate variable from coor in nc 
% latvar - string of the latitude coordinate variable from coor in nc 
%
% If multiple possibilities exist in coor for lat/long coordinates, the
% function throws an error.
%
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
% 
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if debug; 
 nc
end

lonvar=[]; latvar=[];
for i=1:length(coor)
	 name=char(coor{i});
	 if debug; disp(['Variable = ',name]); end
	 if isfield(nc.(name),'units') & ~isempty(strfind(nc.(name).units,'degrees'))
	 	 if ~isempty(strfind(name,'lon')) | ...
		    ~isempty(strfind(name,'Lon')) | ...
		    ~isempty(strfind(name,'LON'))
		    	if ~isempty(lonvar)
			 error(['Two LON candidates in: ',nccellcat(coor)]);
			end
			lonvar=name;
		 end 
	 	 if ~isempty(strfind(name,'lat')) | ...
		    ~isempty(strfind(name,'Lat')) | ...
		    ~isempty(strfind(name,'LAT'))
		    	if ~isempty(lonvar)
			 error(['Two LAT candidates in: ',nccellcat(coor)]);
			end
			latvar=name;
		 end 
	 elseif strcmp(name,'longitude')
	  	disp('Assuming latitude units are degrees_East')
	  	nc.longitude.units='degrees_East';
	  	lonvar=name;
	 elseif strcmp(name,'latitude')
	  	disp('Assuming latitude units are degrees_North')
	  	nc.latitude.units='degrees_North';
 	  	latvar=name;
	 end
end
if isempty(lonvar) | isempty(lonvar)
	 error(['Can''t find lat/lon coordinates in: ',nccellcat(coor)]);
end

if debug; disp(['...Leaving ',mfilename]); end

