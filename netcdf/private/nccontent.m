function [dimc,coor,dims,vars,mdim,mesh,supp]=nccontent(nc)

% NCCONTENT - Classifies nc data into dimensions, coordinates, variables and supporting variables
%
% function [dimc,coor,dims,vars,mdim,mesh,supp]=nccontent(nc)
%
% This function defines the variables in a netCDF structure and their
% purpose in the dataset. It is primarily of use when data is being 
% manipulated and dimensions or data are being removed during analysis.
%
% Input:
% nc - netCDF structure (see help ncstruct for more information).
% 
% Output defining the content of nc:
% dimc - cell array of dimensions used to define coordinates.
% coor - cell array of coordinates used to geolocate variables.
% dims - cell array of dimensions used to define variables.
% vars - cell array of variables defining main data content.
% mdim - cell array of dimensions used to define a mesh and not a variable.
% mesh - cell array of variables defining a mesh but not other data.
% supp - cell array of supporting variables
%
% Example:
% Given the following nc structure:
%
% netCDF structure nc {
%  attributes {global}
%     Conventions: 'CF-1.0'
%         history: '2007-11-12 20:00:39 GMT by mars2netcdf-0.92'
%           title: 'era40.1958'
%  
%  5 variables, dimension=[ ], data=( )
%  
%  latitude [-90 to 90, range 180, step -2.5]
%          data: [73x1 double]
%     dimension: {'latitude'}
%     long_name: 'latitude'
%          type: 'NC_FLOAT'
%         units: 'degrees_north'
%  
%  longitude [0 to 357.5, range 357.5, step 2.5]
%          data: [1x144 double]
%     dimension: {'longitude'}
%     long_name: 'longitude'
%          type: 'NC_FLOAT'
%         units: 'degrees_east'
%  
%  time [01-Jan-1958 00:00:00 to 01-Jan-1959 00:00:00, timestep 06:00:00.000]
%          data: [1461x1 double]
%     dimension: {'time'}
%     long_name: 'time'
%          type: 'NC_FLOAT'
%         units: 'hours since 1900-01-01 00:00:0.0'
%        weight: [1x1461 double]
%  
%  cell_area (8.2564e-07 to 0.00015166, range 0.00015083, median 0.0001071)
%     coordinates: 'latitude longitude'
%            data: [73x144 double]
%       dimension: {'latitude'  'longitude'}
%       long_name: 'fractional area on the sphere'
%            type: 'NC_DOUBLE'
%  
%  p2t (193.1944 to 321.9842, range 128.7898, median 282.9082)
%     coordinates: 'latitude longitude'
%            data: [73x144x1461 double]
%       dimension: {'latitude'  'longitude'  'time'}
%       fillvalue: -99999
%       long_name: '2 metre temperature'
%            type: 'NC_FLOAT'
%           units: 'K'
%  
% }
%
% The command "[dimc,coor,dims,vars,mdim,mesh]=nccontent(nc)" yields:
%
% Checking the content of the netCDF structure
%   Coordinates set for: cell_area
%   Coordinates set for: p2t
% Coordinate dimensions: latitude longitude
%  Coordinate variables: latitude longitude
%   Variable dimensions: time
%        Data variables: p2t
%  Supporting variables: cell_area
%  
% dimc = 
%     'latitude'    'longitude'
% coor = 
%     'latitude'    'longitude'
% dims = 
%     'time'
% vars = 
%     'p2t'
% mdim = 
%      {}
% mesh = 
%      {}
%
% Coordinate variables are defined by the attribute 'coordinates' such as:
% nc.p2t.coordinates='latitude longitude'
%
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
% 
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if ~isstruct(nc); error('nc input is not a structure'); end

if debug; disp('Checking the content of the netCDF structure'); end

% first, pass through each variable of the netCDF structure to find dimensions and coordinates 
coor={};
dims={};
vars={};
mesh={};
supp={};

[nc,variablenames,numbervariables]=ncsort(nc);
for m = 1:numbervariables

	 var=char(variablenames(m));

	 % find coordinates
	 varcoords=ncvarcoords(nc,var);
         
	 if length(varcoords)>0
		coor=ncunion(coor,varcoords);
	 end

	 % find mesh variables
	 if isfield(nc.(var),'units') & ~any(strcmp(var,nc.(var).dimension))
 	    if strcmp(char(nc.(var).units),'degrees_east') | ...
	       strcmp(char(nc.(var).units),'degrees_west') | ...
	       strcmp(char(nc.(var).units),'degrees_north') | ...
	       strcmp(char(nc.(var).units),'degrees_south') 

                varcoords={var};
		mesh=ncunion(mesh,varcoords);

	    end
	 end

	 % find support grid information
	 if isfield(nc.(var),'grid_mapping_name')
                supp={var};
	 end

	 % find data variables 
	 if ~any(strcmp(var,nc.(var).dimension))
                varcoords={var};
		vars=ncunion(vars,varcoords);
	 end

end


% remove dimension variables from other sets
vars=ncsetdiff(vars,mesh);
vars=ncsetdiff(vars,coor);
mesh=ncsetdiff(mesh,coor);

% define supporting material, including search for grid information
supporting={'cell_area',...
	    'time_bounds',...
	    'time_string',...
	    'mask',...
            'turn',...
            'hmask',...
            'tmask',...
            'rotang',...
            'tarea',...
            'uarea',...
            'dxt',...
            'dyt',...
            'dxu',...
            'dyu',...
            'blkmask',...
            'ANGLE',...
            'ANGLET',...
            'DXT',...
            'DXU',...
            'DYT',...
            'DYU',...
            'HT',...
            'HTE',...
            'HTN',...
            'HU',...
            'HUS',...
            'HUW',...
            'KMT',...
            'KMU',...
            'TAREA',...
            'UAREA'};

% remove supporting variables from the variable list to supp
supp=ncunion(supp,intersect(supporting,vars));
vars=ncsetdiff(vars,supp);

% Now combine these unique sets to define the dataset
% Get dimensions that define variables
dims={};
for i=1:length(vars)
 dims=ncunion(dims,nc.(char(vars{i})).dimension);
end

% check coordinate dimensions are a subset or equal set to dims
dimc={};
coordcrash=[];
for i=1:length(coor)
 if isfield(nc,char(coor{i}))      
  dimc=ncunion(dimc,nc.(char(coor{i})).dimension);
 else
  coordcrash=[coordcrash,char(coor{i}),' '];
 end
end
if ~isempty(coordcrash)
 error(['Please also import coordinates: ',coordcrash])
end


% remove dimc dimensions from dims (SOMETHING MAY BE WRONG HERE BY 
% REMOVING THIS LINE - POSSIBLE BUG)
% dims=ncsetdiff(dims,dimc);

%  Get mesh dimensions that define the mesh and are not in dims
mdim={};
for i=1:length(mesh)
	mdim=ncunion(mdim,nc.(char(mesh{i})).dimension);
end
mdim=ncsetdiff(mdim,ncunion(dimc,dims));

% output information and checks
if ~isempty(dimc) & nargout<1 
	disp(['Coordinate dimensions:',nccellcat(dimc)])
end

if ~isempty(coor) & nargout<1
	disp([' Coordinate variables:',nccellcat(coor)])
end

if isempty(dims) 
	if debug; disp('WARNING: No variable dimensions detected'); end
elseif nargout<1
	disp(['  Variable dimensions:',nccellcat(dims)])
end

if isempty(vars) 
	if debug; disp('WARNING: No data variables have been detected in this dataset'); end
elseif nargout<1
	disp(['       Data variables:',nccellcat(vars)])
end

if ~isempty(supp) & nargout<1
	if debug; disp([' Supporting variables:',nccellcat(supp)]); end
end

if ~isempty(mesh) & nargout<1
	disp(['      Mesh dimensions:',nccellcat(mdim)])
	disp(['       Mesh variables:',nccellcat(mesh)])
end

if debug; disp(' '); end

if debug; disp(['...Leaving ',mfilename]); end

