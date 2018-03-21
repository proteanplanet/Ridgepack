function [gridtype]=ridgepack_gridtype(nc)

% ridgepack_gridtype - Determines if a grid is regional or global
%
% function [gridtype]=ridgepack_gridtype(nc)
%
% This function determines if this is a global or regional 
% model grid based on a few simple assumptions about the 
% extent of latitude and longitude.
%
% INPUT:
%
% nc - nc structure containing the 2D fields 'latitude' and 'longitude'
%      
%
% OUTPUT:
%
% gridtype - a text string as either 'regional' or 'global'
%
% Note that nc.attributes.gridtype can be manually set before calling
% this function to explicitly set the type, rather than relying on automation.
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

% Check that nc is a structure
if not(isstruct(nc));
   error([inputname(1),' is not a structure']);
end

% Generate 2D lats and longs in case this has not already been done.
nc=ridgepack_sph2gen(nc);

% determine if the grid is global or regional
gridtype='nonsense';

if isfield(nc.attributes,'gridtype')
 gridtype=nc.attributes.gridtype;
else
 if isfield(nc,'longitude') & isfield(nc,'latitude') & ...
   min(size(nc.longitude.data))>1 & min(size(nc.latitude.data))>1 & ...
   size(nc.longitude.data)==size(nc.latitude.data)
   if  (range(nc.longitude.data(:))>350  & range(nc.latitude.data(:))>120)
        gridtype='global';
   elseif  (range(nc.longitude.data(:))<350 | range(nc.latitude.data(:))<90)
        gridtype='regional';
   else
    while ~any(strcmp(gridtype,{'regional','global'}))
      gridtype=input('Enter ''regional'' or ''global'' dataset:','s');
    end
   end
 elseif size(nc.longitude.data)~=size(nc.latitude.data)
   error('latitude is a different matrix size to longitude');
 elseif isfield(nc,'longitude') & isfield(nc,'latitude') 
   error('latitude and/or longitude are not 2D fields');
 else
   error('Unable to assign gridtype');
 end
end

if ~any(strcmp(gridtype,{'regional','global'}))
   error('Grid type neither regional nor global')
else
   if debug; disp(['Dataset is ',gridtype]); end
end

if debug; disp(['...Leaving ',mfilename]); end

