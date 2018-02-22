function nc=ridgepack_geodetictoxy_icetide(nc)

% ridgepack_geodetictoxy_ICETIDE - Convert lat-long to X and Y coordinates on the ice-tide grid
%
% function nc=ridgepack_geodetictoxy_icetide(nc)
%
% This function converts latitudes and longitudes to x and y coordinates
% on the Hibler and Roberts ice-tide grid.
% 
% INPUT:
%
% nc - nc structure containing lats and longs from the model
%
%
% OUTPUT:
%
% nc - nc structure containing lats and longs from the model, and model positions
%
% Ridgepack Version 1.0
% Andrew Roberts, Naval Postgraduate School, March 2018 (afrobert@nps.edu)
%

global debug;
if debug; disp(['Entering ',mfilename,'...']); end

if ~isfield(nc,'latitude')
 error('latitude is missing from the nc structure')
elseif ~isfield(nc,'longitude')
 error('longitude is missing from the nc structure')
end

r=8.0*(90.0-nc.latitude.data);
nc.x.data=r.*cos((pi./180).*(nc.longitude.data-40.0))+180.0;
nc.y.data=r.*sin((pi./180).*(nc.longitude.data-40.0))+133.0;

nc.x.dimension=nc.latitude.dimension;
nc.y.dimension=nc.latitude.dimension;
nc.x.long_name='Model x-position';
nc.y.long_name='Model y-position';

nc=ridgepack_struct(nc);

if debug; disp(['...Leaving ',mfilename]); end

