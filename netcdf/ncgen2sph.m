function [nc]=ncgen2sph(nc)

% NCGEN2SPH - Convert nc structure from general coordinates to lat-lon dimensions
%
% function [nc]=ncgen2sph(nc)
%
% This function works in reverse to the function ncsph2gen by taking an nc structure
% in standard generalized format and putting into spherical coordinate structure
% if it is indeed possible.  To understand the structures. see the explanation in
% the help page for ncsph2gen.
%
% Input:
% nc - netcdf structure (see ncstruct for more information) with data 
%      placed in a standardarised generalised coodinates on a sphere.  
%
% Output:
% nc - netcdf structure (see ncstruct for more information) with data in 
%      spherical coordinates.
%
% Written by Andrew Roberts 2012
% Department of Oceanography, Naval Postgraduate School
%
% $Id$

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if debug; disp('Generating lat-lon grid from generalized grid'); end

% Check that nc is a structure
if not(isstruct(nc));
 error([inputname(1),' is not a structure']);
end

% Check that nc is in general coordinates
if ~isfield(nc,'x') | ~isfield(nc,'y')
 disp('not in general coordinate format for ncbox')
 error('x and/or y are not fields')
elseif ~strcmp('x',nc.x.dimension) | ~strcmp('y',nc.y.dimension)
 disp('not in general coordinate format for ncbox')
 error('The dimensions of x and/or y are incorrect')
elseif ~(length(nc.longitude.dimension)>1 & ~strcmp('longitude',nc.longitude.dimension) & ...
         length(nc.latitude.dimension)>1 & ~strcmp('latitude',nc.latitude.dimension))
 disp('not in general coordinate format for ncbox')
 error('This nc structure does not have 2D lat and long arrays')
end

% Check that general coordinates can be put in spherical coordinates
if nc.longitude.data(1,:)~=nc.longitude.data(2,:)
 disp('The longitude values in row 1 differ from row 2')
 error('The generalised coordinates cannot be converted to spherical')
elseif nc.latitude.data(:,1)~=nc.latitude.data(:,2)
 disp('The latitude values in column 1 differ from column 2')
 error('The generalised coordinates cannot be converted to spherical')
end

if isfield(nc,'longitude_corner') & ...
   nc.longitude_corner.data(1,:)~=nc.longitude_corner.data(2,:)
 disp('The longitude_corner values in row 1 differ from row 2')
 error('The generalised coordinates cannot be converted to spherical')
elseif isfield(nc,'latitude_corner') & ...
   nc.latitude_corner.data(:,1)~=nc.latitude_corner.data(:,2)
 disp('The latitude_corner values in column 1 differ from column 2')
 error('The generalised coordinates cannot be converted to spherical')
end
 
% Change the dimension in other variables
[nc,variablenames,numbervariables]=ncsort(nc);
for m = 1:numbervariables
 name=char(variablenames(m));
 for i=1:length(nc.(name).dimension)
  if strcmp('x',char(nc.(name).dimension{i}))
      nc.(name).dimension{i}='longitude';
  end
  if strcmp('y',char(nc.(name).dimension{i}))
      nc.(name).dimension{i}='latitude';
  end
  if strcmp('x_corner',char(nc.(name).dimension{i}))
      nc.(name).dimension{i}='longitude_corner';
  end
  if strcmp('y_corner',char(nc.(name).dimension{i}))
      nc.(name).dimension{i}='latitude_corner';
  end
 end
end

% Remove x and y dimensions
nc=rmfield(nc,'x');
nc=rmfield(nc,'y');

% Remove x_corner and y_corner dimensions if they exist
if isfield(nc,'x_corner'); nc=rmfield(nc,'x_corner'); end;
if isfield(nc,'y_corner'); nc=rmfield(nc,'y_corner'); end;

% Make latitudes and longitudes the dimension variables
lon=nc.longitude.data(1,:); nc.longitude.data=lon;
nc.longitude.dimension={'longitude'};
lat=nc.latitude.data(:,1); nc.latitude.data=lat;
nc.latitude.dimension={'latitude'};

% Make latitudes_corner and longitudes_corner the dimension variables if they exist
if isfield(nc,'longitude_corner'); 
 lon=nc.longitude_corner.data(1,:); nc.longitude_corner.data=lon;
 nc.longitude_corner.dimension={'longitude_corner'};
end
if isfield(nc,'latitude_corner'); 
 lat=nc.latitude_corner.data(:,1); nc.latitude_corner.data=lat;
 nc.latitude_corner.dimension={'latitude_corner'};
end

% run check through the data
[nc,out]=ncstruct(nc);

if debug; disp(['...Leaving ',mfilename]); end

